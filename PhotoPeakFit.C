// ROOT macro for RadWare/GF3-style single photopeak fitting.

#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TH1.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TObject.h"
#include "TROOT.h"
#include "TString.h"
#include "TVirtualPad.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>

// Default values copied from RadWare/GF3 where applicable:
// R = 10, BETA = W0 / 2, STEP = 0.25,
// W0 = sqrt(9 + 0.004 * peak0).
// Histogram-dependent defaults are set inside photopeakfit():
// bg0, bg1, bg2 from the fit endpoints; Centroid from peak0;
// Height from peak0 minus the estimated linear background.
enum PhotoPeakPar {
  kA = 0,     // bg0: constant term of the quadratic background.
  kB = 1,     // bg1: linear term of the quadratic background.
  kC = 2,     // bg2: quadratic term of the quadratic background.
  kR = 3,     // R: skew-tail fraction in percent.
  kBeta = 4,  // BETA: skew-tail decay constant.
  kStep = 5,  // STEP: smoothed step relative height.
  kP = 6,     // Centroid: peak centroid position.
  kW = 7,     // FWHM: full width at half maximum.
  kH = 8,     // Height: fitted peak height.
  kNPars = 9  // Number of fit parameters.
};

enum PhotoPeakFitMode {
  kPhotoPeakAuto = 0,
  kPhotoPeakHighStat = 1,
  kPhotoPeakLowStat = 2
};

struct PhotoPeakFitCandidate {
  const char *name;
  bool useTail;
  bool useStep;
  bool useQuadBg;
};

struct PhotoPeakFitTrial {
  TF1 *func;
  TMatrixDSym cov;
  PhotoPeakFitCandidate candidate;
  int status;
  int ndf;
  int freePars;
  double chi2;
  double reducedChi2;
};

struct PhotoPeakFitConfig {
  int mode;
  bool hasInit[kNPars];
  double init[kNPars];
  bool hasLimit[kNPars];
  double limitLow[kNPars];
  double limitHigh[kNPars];
  bool hasFix[kNPars];
  double fix[kNPars];

  PhotoPeakFitConfig()
  {
    mode = kPhotoPeakAuto;
    for (int i = 0; i < kNPars; ++i) {
      hasInit[i] = false;
      init[i] = 0.0;
      hasLimit[i] = false;
      limitLow[i] = 0.0;
      limitHigh[i] = 0.0;
      hasFix[i] = false;
      fix[i] = 0.0;
    }
  }
};

double gPhotoPeakFitLow = 0.0;
double gPhotoPeakFitHigh = 0.0;

int PhotoPeakFreeParameterCount(const PhotoPeakFitCandidate &candidate);

// ============== PhotoPeakTrim ==============
// Purpose: Remove leading and trailing whitespace from a string.
// Inputs: String value.
// Outputs: Trimmed string.
std::string PhotoPeakTrim(const std::string &text)
{
  size_t first = 0;
  while (first < text.size() &&
         std::isspace(static_cast<unsigned char>(text[first]))) {
    ++first;
  }
  size_t last = text.size();
  while (last > first &&
         std::isspace(static_cast<unsigned char>(text[last - 1]))) {
    --last;
  }

  return text.substr(first, last - first);
}

// ============== PhotoPeakLower ==============
// Purpose: Convert a string to lowercase for keyword matching.
// Inputs: String value.
// Outputs: Lowercase string.
std::string PhotoPeakLower(const std::string &text)
{
  std::string out = text;
  for (size_t i = 0; i < out.size(); ++i) {
    out[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(out[i])));
  }

  return out;
}

// ============== PhotoPeakUpper ==============
// Purpose: Convert a string to uppercase for parameter matching.
// Inputs: String value.
// Outputs: Uppercase string.
std::string PhotoPeakUpper(const std::string &text)
{
  std::string out = text;
  for (size_t i = 0; i < out.size(); ++i) {
    out[i] = static_cast<char>(std::toupper(static_cast<unsigned char>(out[i])));
  }

  return out;
}

// ============== PhotoPeakParName ==============
// Purpose: Convert a parameter index to the public parameter name.
// Inputs: Parameter index.
// Outputs: Static parameter name string.
const char *PhotoPeakParName(int ipar)
{
  const char *names[kNPars] = {
    "A", "B", "C", "R", "BETA", "STEP", "P", "W", "H"
  };
  if (ipar < 0 || ipar >= kNPars) {
    return "UNKNOWN";
  }

  return names[ipar];
}

// ============== PhotoPeakParIndex ==============
// Purpose: Convert a public parameter name to its parameter index.
// Inputs: Parameter name string.
// Outputs: Parameter index, or -1 when unknown.
int PhotoPeakParIndex(const std::string &name)
{
  const std::string key = PhotoPeakUpper(PhotoPeakTrim(name));
  for (int i = 0; i < kNPars; ++i) {
    if (key == PhotoPeakParName(i)) {
      return i;
    }
  }

  return -1;
}

// ============== PhotoPeakEval ==============
// Purpose: Evaluate the full RadWare/GF3-style photopeak fit function.
// Inputs: ROOT x array and parameter array.
// Outputs: Total fit value at x.
double PhotoPeakEval(double *x, double *par)
{
  const double xx = x[0];
  const double a = par[kA];
  const double b = par[kB];
  const double c = par[kC];
  const double rPercent = par[kR];
  const double beta = std::max(par[kBeta], 1.0e-12);
  const double step = par[kStep];
  const double p = par[kP];
  const double wFwhm = std::max(par[kW], 1.0e-12);
  const double h = par[kH];

  const double sigma = wFwhm / 2.35482;
  const double w = (xx - p) / (sigma * TMath::Sqrt2());
  const double y = wFwhm / (beta * 3.33021838);
  const double erfcY = std::max(TMath::Erfc(y), 1.0e-300);
  const double gaussian = std::exp(-w * w);
  const double stepShape = TMath::Erfc(w);
  double skew = 0.0;

  const double tailArg = (xx - p) / beta;
  if (std::fabs(tailArg) <= 700.0) {
    skew = std::exp(tailArg) * TMath::Erfc(w + y) / erfcY;
  }

  const double r = rPercent / 100.0;
  const double xCentered = xx - 0.5 * (gPhotoPeakFitLow + gPhotoPeakFitHigh);
  const double quadBg = a + b * xCentered + c * xCentered * xCentered;
  const double stepBg = h * step * stepShape / 200.0;
  const double photopeak = h * ((1.0 - r) * gaussian + r * skew);

  return photopeak + quadBg + stepBg;
}

// ============== PhotoPeakBgEval ==============
// Purpose: Evaluate the background-only part of the fit.
// Inputs: ROOT x array and parameter array.
// Outputs: Quadratic plus smoothed step background value at x.
double PhotoPeakBgEval(double *x, double *par)
{
  const double xx = x[0];
  const double sigma = std::max(par[kW], 1.0e-12) / 2.35482;
  const double w = (xx - par[kP]) / (sigma * TMath::Sqrt2());
  const double xCentered = xx - 0.5 * (gPhotoPeakFitLow + gPhotoPeakFitHigh);
  const double quadBg = par[kA] + par[kB] * xCentered + par[kC] * xCentered * xCentered;
  const double stepBg = par[kH] * par[kStep] * TMath::Erfc(w) / 200.0;

  return quadBg + stepBg;
}

// ============== PhotoPeakBinContentAtX ==============
// Purpose: Read histogram content at the bin containing an x-value.
// Inputs: Histogram pointer and x-value.
// Outputs: Bin content.
double PhotoPeakBinContentAtX(TH1 *hist, double x)
{
  const int bin = hist->GetXaxis()->FindFixBin(x);
  const int clampedBin = std::max(1, std::min(hist->GetNbinsX(), bin));

  return hist->GetBinContent(clampedBin);
}

// ============== PhotoPeakMaxInRange ==============
// Purpose: Find the largest bin content inside an x range.
// Inputs: Histogram pointer and x-value limits.
// Outputs: Maximum bin content in the range.
double PhotoPeakMaxInRange(TH1 *hist, double fitLow, double fitHigh)
{
  int lowBin = hist->GetXaxis()->FindFixBin(fitLow);
  int highBin = hist->GetXaxis()->FindFixBin(fitHigh);
  lowBin = std::max(1, std::min(hist->GetNbinsX(), lowBin));
  highBin = std::max(1, std::min(hist->GetNbinsX(), highBin));
  if (lowBin > highBin) {
    std::swap(lowBin, highBin);
  }

  double maxContent = hist->GetBinContent(lowBin);
  for (int bin = lowBin + 1; bin <= highBin; ++bin) {
    maxContent = std::max(maxContent, hist->GetBinContent(bin));
  }

  return maxContent;
}

// ============== PhotoPeakBinCountInRange ==============
// Purpose: Count histogram bins touched by an x range.
// Inputs: Histogram pointer and x-value limits.
// Outputs: Inclusive bin count inside the range.
int PhotoPeakBinCountInRange(TH1 *hist, double fitLow, double fitHigh)
{
  int lowBin = hist->GetXaxis()->FindFixBin(fitLow);
  int highBin = hist->GetXaxis()->FindFixBin(fitHigh);
  lowBin = std::max(1, std::min(hist->GetNbinsX(), lowBin));
  highBin = std::max(1, std::min(hist->GetNbinsX(), highBin));
  if (lowBin > highBin) {
    std::swap(lowBin, highBin);
  }

  return highBin - lowBin + 1;
}

// ============== PhotoPeakModeName ==============
// Purpose: Convert a user mode value to a printable name.
// Inputs: Fit mode integer.
// Outputs: Static mode name string.
const char *PhotoPeakModeName(int mode)
{
  if (mode == kPhotoPeakHighStat) {
    return "highstat";
  }
  if (mode == kPhotoPeakLowStat) {
    return "lowstat";
  }

  return "auto";
}

// ============== PhotoPeakParseMode ==============
// Purpose: Convert a user mode string to a fitting mode.
// Inputs: Mode string and output mode reference.
// Outputs: True when the mode string is recognized.
bool PhotoPeakParseMode(const std::string &text, int &mode)
{
  const std::string key = PhotoPeakLower(PhotoPeakTrim(text));
  if (key == "auto" || key == "kphotopeakauto") {
    mode = kPhotoPeakAuto;
    return true;
  }
  if (key == "highstat" || key == "high_stat" ||
      key == "kphotopeakhighstat") {
    mode = kPhotoPeakHighStat;
    return true;
  }
  if (key == "lowstat" || key == "low_stat" ||
      key == "kphotopeaklowstat") {
    mode = kPhotoPeakLowStat;
    return true;
  }

  return false;
}

// ============== PhotoPeakCandidateUsesPar ==============
// Purpose: Check whether a candidate function uses a parameter.
// Inputs: Fit candidate and parameter index.
// Outputs: True when the parameter is active for the candidate.
bool PhotoPeakCandidateUsesPar(const PhotoPeakFitCandidate &candidate, int ipar)
{
  if (ipar == kA || ipar == kB || ipar == kP ||
      ipar == kW || ipar == kH) {
    return true;
  }
  if (ipar == kC) {
    return candidate.useQuadBg;
  }
  if (ipar == kR || ipar == kBeta) {
    return candidate.useTail;
  }
  if (ipar == kStep) {
    return candidate.useStep;
  }

  return false;
}

// ============== PhotoPeakConfiguredFreeParameterCount ==============
// Purpose: Count candidate free parameters after user fixed parameters.
// Inputs: Fit candidate and configuration.
// Outputs: Number of free parameters.
int PhotoPeakConfiguredFreeParameterCount(const PhotoPeakFitCandidate &candidate,
                                          const PhotoPeakFitConfig &config)
{
  int freePars = PhotoPeakFreeParameterCount(candidate);
  for (int ipar = 0; ipar < kNPars; ++ipar) {
    if (config.hasFix[ipar] && PhotoPeakCandidateUsesPar(candidate, ipar)) {
      --freePars;
    }
  }

  return std::max(0, freePars);
}

// ============== PhotoPeakReadConfig ==============
// Purpose: Read a PhotoPeakFit text configuration file.
// Inputs: Configuration file path and output configuration reference.
// Outputs: Zero on success, non-zero on file open error.
int PhotoPeakReadConfig(const char *configFile, PhotoPeakFitConfig &config)
{
  if (!configFile || !configFile[0]) {
    std::printf("photopeakfit ERROR: empty config file path.\n");
    return 1;
  }

  std::ifstream input(configFile);
  if (!input) {
    std::printf("photopeakfit ERROR: cannot open config file '%s'.\n", configFile);
    return 2;
  }

  std::string line;
  int lineNumber = 0;
  while (std::getline(input, line)) {
    ++lineNumber;
    const size_t commentPos = line.find('#');
    if (commentPos != std::string::npos) {
      line = line.substr(0, commentPos);
    }
    for (size_t i = 0; i < line.size(); ++i) {
      if (line[i] == '=') {
        line[i] = ' ';
      }
    }
    line = PhotoPeakTrim(line);
    if (line.empty()) {
      continue;
    }

    std::istringstream words(line);
    std::string keyword;
    words >> keyword;
    keyword = PhotoPeakLower(keyword);

    if (keyword == "mode") {
      std::string modeText;
      words >> modeText;
      int parsedMode = kPhotoPeakAuto;
      if (!PhotoPeakParseMode(modeText, parsedMode)) {
        std::printf("photopeakfit WARNING: %s:%d unknown mode '%s'; using auto.\n",
                    configFile, lineNumber, modeText.c_str());
        parsedMode = kPhotoPeakAuto;
      }
      config.mode = parsedMode;
      continue;
    }

    if (keyword == "init" || keyword == "fix" || keyword == "limit") {
      std::string parName;
      words >> parName;
      const int ipar = PhotoPeakParIndex(parName);
      if (ipar < 0) {
        std::printf("photopeakfit WARNING: %s:%d unknown parameter '%s'; ignoring line.\n",
                    configFile, lineNumber, parName.c_str());
        continue;
      }

      if (keyword == "init" || keyword == "fix") {
        double value = 0.0;
        if (!(words >> value)) {
          std::printf("photopeakfit WARNING: %s:%d missing value; ignoring line.\n",
                      configFile, lineNumber);
          continue;
        }
        if (keyword == "init") {
          config.hasInit[ipar] = true;
          config.init[ipar] = value;
        } else {
          config.hasFix[ipar] = true;
          config.fix[ipar] = value;
        }
        continue;
      }

      double low = 0.0;
      double high = 0.0;
      if (!(words >> low >> high)) {
        std::printf("photopeakfit WARNING: %s:%d missing limit values; ignoring line.\n",
                    configFile, lineNumber);
        continue;
      }
      if (low > high) {
        std::swap(low, high);
      }
      config.hasLimit[ipar] = true;
      config.limitLow[ipar] = low;
      config.limitHigh[ipar] = high;
      continue;
    }

    std::printf("photopeakfit WARNING: %s:%d unknown keyword '%s'; ignoring line.\n",
                configFile, lineNumber, keyword.c_str());
  }

  return 0;
}

// ============== PhotoPeakFreeParameterCount ==============
// Purpose: Count free parameters for a candidate fitting function.
// Inputs: Fit candidate description.
// Outputs: Number of free parameters.
int PhotoPeakFreeParameterCount(const PhotoPeakFitCandidate &candidate)
{
  int freePars = 5;
  if (candidate.useTail) {
    freePars += 2;
  }
  if (candidate.useStep) {
    freePars += 1;
  }
  if (candidate.useQuadBg) {
    freePars += 1;
  }

  return freePars;
}

// ============== PhotoPeakConfigureFunction ==============
// Purpose: Set initial values, limits, and fixed parameters for one fit.
// Inputs: Function, candidate, initial values, and fit limits.
// Outputs: Configured ROOT function.
void PhotoPeakConfigureFunction(TF1 *func, const PhotoPeakFitCandidate &candidate,
                                double a0, double b0, double c0,
                                double r0, double beta0, double step0,
                                double peak0, double w0, double h0,
                                double fitLow, double fitHigh,
                                double range, double hUpper,
                                const PhotoPeakFitConfig &config)
{
  func->SetParNames("A", "B", "C", "R", "BETA", "STEP", "P", "W", "H");
  func->SetParameters(a0, b0, c0, r0, beta0, step0, peak0, w0, h0);
  func->SetParLimits(kR, 0.0, 100.0);
  func->SetParLimits(kBeta, 1.0e-6, 10.0 * range);
  func->SetParLimits(kStep, 0.0, 100.0);
  func->SetParLimits(kP, fitLow, fitHigh);
  func->SetParLimits(kW, 1.0e-6, range);
  func->SetParLimits(kH, 0.0, hUpper);
  if (!candidate.useTail) {
    func->FixParameter(kR, 0.0);
    func->FixParameter(kBeta, beta0);
  }
  if (!candidate.useStep) {
    func->FixParameter(kStep, 0.0);
  }
  if (!candidate.useQuadBg) {
    func->FixParameter(kC, 0.0);
  }
  for (int ipar = 0; ipar < kNPars; ++ipar) {
    if (!PhotoPeakCandidateUsesPar(candidate, ipar)) {
      continue;
    }
    if (config.hasInit[ipar]) {
      func->SetParameter(ipar, config.init[ipar]);
    }
    if (config.hasLimit[ipar]) {
      func->SetParLimits(ipar, config.limitLow[ipar], config.limitHigh[ipar]);
    }
    if (config.hasFix[ipar]) {
      func->FixParameter(ipar, config.fix[ipar]);
    }
  }
  func->SetNpx(2000);
}

// ============== PhotoPeakWarnInactiveConfig ==============
// Purpose: Warn when config lines target parameters inactive in the final model.
// Inputs: Final fit candidate and configuration.
// Outputs: Warning messages printed to stdout.
void PhotoPeakWarnInactiveConfig(const PhotoPeakFitCandidate &candidate,
                                 const PhotoPeakFitConfig &config)
{
  for (int ipar = 0; ipar < kNPars; ++ipar) {
    if (PhotoPeakCandidateUsesPar(candidate, ipar)) {
      continue;
    }
    if (config.hasInit[ipar]) {
      std::printf("photopeakfit WARNING: init %s is inactive in %s; ignoring.\n",
                  PhotoPeakParName(ipar), candidate.name);
    }
    if (config.hasLimit[ipar]) {
      std::printf("photopeakfit WARNING: limit %s is inactive in %s; ignoring.\n",
                  PhotoPeakParName(ipar), candidate.name);
    }
    if (config.hasFix[ipar]) {
      std::printf("photopeakfit WARNING: fix %s is inactive in %s; ignoring.\n",
                  PhotoPeakParName(ipar), candidate.name);
    }
  }
}

// ============== PhotoPeakTrialIsBetter ==============
// Purpose: Decide whether one auto-mode trial is better than another.
// Inputs: New trial pointer and current best trial pointer.
// Outputs: True when the new trial should become the best fit.
bool PhotoPeakTrialIsBetter(const PhotoPeakFitTrial *trial, const PhotoPeakFitTrial *best)
{
  if (!trial || trial->ndf <= 0 || !std::isfinite(trial->reducedChi2)) {
    return false;
  }
  if (!best || best->ndf <= 0 || !std::isfinite(best->reducedChi2)) {
    return true;
  }
  if (trial->status == 0 && best->status != 0) {
    return true;
  }
  if (trial->status != 0 && best->status == 0) {
    return false;
  }

  return trial->reducedChi2 < best->reducedChi2;
}

// ============== PhotoPeakArea ==============
// Purpose: Calculate RadWare/GF3 photopeak area from the photopeak terms only.
// Inputs: Fit parameter array and histogram bin width.
// Outputs: Photopeak area in histogram-bin counts.
double PhotoPeakArea(const double *par, double binWidth)
{
  const double r = par[kR] / 100.0;
  const double beta = std::max(par[kBeta], 1.0e-12);
  const double wFwhm = std::max(par[kW], 1.0e-12);
  const double h = par[kH];
  const double y = wFwhm / (beta * 3.33021838);
  const double erfcY = TMath::Erfc(y);
  const double d = erfcY > 0.0 ? std::exp(-y * y) / erfcY : 0.0;

  const double continuousArea = h * (r * 2.0 * beta * d + (1.0 - r) * wFwhm * 1.06446705);
  return continuousArea / std::max(binWidth, 1.0e-12);
}

// ============== PhotoPeakAreaUncertainty ==============
// Purpose: Propagate fit covariance to the photopeak area uncertainty.
// Inputs: Fit parameter array, covariance matrix, and histogram bin width.
// Outputs: One-sigma area uncertainty in histogram-bin counts.
double PhotoPeakAreaUncertainty(const double *par, const TMatrixDSym &cov, double binWidth)
{
  const double invBinWidth = 1.0 / std::max(binWidth, 1.0e-12);
  const double r = par[kR] / 100.0;
  const double beta = std::max(par[kBeta], 1.0e-12);
  const double wFwhm = std::max(par[kW], 1.0e-12);
  const double h = par[kH];
  const double y = wFwhm / (beta * 3.33021838);
  const double erfcY = TMath::Erfc(y);
  const double d = erfcY > 0.0 ? std::exp(-y * y) / erfcY : 0.0;
  const double dAdH = r * 2.0 * beta * d + (1.0 - r) * wFwhm * 1.06446705;
  const double dAdR = 0.01 * h * (2.0 * beta * d - wFwhm * 1.06446705);
  const double dAdBeta = h * r * 2.0 * d *
    (1.0 + 2.0 * y * y - d * 1.12837917 * y);
  const double dAdW = h *
    (r * 2.0 * 0.600561216 * d * (d / 1.77245385 - y) +
     (1.0 - r) * 1.06446705);
  double grad[kNPars] = {0.0};
  grad[kR] = dAdR * invBinWidth;
  grad[kBeta] = dAdBeta * invBinWidth;
  grad[kW] = dAdW * invBinWidth;
  grad[kH] = dAdH * invBinWidth;

  double variance = 0.0;
  for (int i = 0; i < kNPars; ++i) {
    for (int j = 0; j < kNPars; ++j) {
      variance += grad[i] * cov(i, j) * grad[j];
    }
  }

  return variance > 0.0 ? std::sqrt(variance) : 0.0;
}

// ============== PhotoPeakFitRun ==============
// Purpose: Fit one histogram photopeak with a selected RadWare/GF3-style function.
// Inputs: Histogram pointer, lower x range, upper x range, initial peak x, and config.
// Outputs: Fit status; prints fit results and draws total/background functions.
int PhotoPeakFitRun(TH1 *hist, double fitLow, double fitHigh, double peak0,
                    const PhotoPeakFitConfig &config)
{
  if (!hist) {
    std::printf("photopeakfit ERROR: null histogram pointer.\n");
    return 1;
  }
  if (fitLow == fitHigh) {
    std::printf("photopeakfit ERROR: fitLow and fitHigh are equal.\n");
    return 2;
  }
  if (fitLow > fitHigh) {
    std::swap(fitLow, fitHigh);
  }
  int mode = config.mode;
  if (mode != kPhotoPeakAuto &&
      mode != kPhotoPeakHighStat &&
      mode != kPhotoPeakLowStat) {
    std::printf("photopeakfit WARNING: unknown mode %d; using auto.\n", mode);
    mode = kPhotoPeakAuto;
  }

  gPhotoPeakFitLow = fitLow;
  gPhotoPeakFitHigh = fitHigh;

  const double yLow = PhotoPeakBinContentAtX(hist, fitLow);
  const double yHigh = PhotoPeakBinContentAtX(hist, fitHigh);
  const double yPeak = PhotoPeakBinContentAtX(hist, peak0);
  const double xMid = 0.5 * (fitLow + fitHigh);
  const double range = fitHigh - fitLow;
  const double a0 = 0.5 * (yLow + yHigh);
  const double b0 = (yHigh - yLow) / range;
  const double c0 = 0.0;
  const double r0 = 10.0;
  const double w0 = std::sqrt(std::max(9.0 + 0.004 * peak0, 1.0e-12));
  const double beta0 = 0.5 * w0;
  const double step0 = 0.25;
  const double linearBgAtPeak = a0 + b0 * (peak0 - xMid);
  const double maxInRange = PhotoPeakMaxInRange(hist, fitLow, fitHigh);
  const double h0 = std::max(yPeak - linearBgAtPeak, std::max(maxInRange, 1.0));
  const double hUpper = std::max(10.0 * maxInRange, h0 * 10.0);
  const int nFitBins = PhotoPeakBinCountInRange(hist, fitLow, fitHigh);
  const PhotoPeakFitCandidate allCandidates[] = {
    {"gaussian_linearBg", false, false, false},
    {"gaussian_linearBg_tail", true, false, false},
    {"gaussian_linearBg_step", false, true, false},
    {"gaussian_linearBg_quadBg", false, false, true},
    {"gaussian_linearBg_tail_step", true, true, false},
    {"gaussian_linearBg_tail_quadBg", true, false, true},
    {"gaussian_linearBg_step_quadBg", false, true, true},
    {"gaussian_linearBg_tail_step_quadBg", true, true, true}
  };
  PhotoPeakFitTrial trials[8] = {
    {0, TMatrixDSym(kNPars), allCandidates[0], 1, 0, 0, 0.0, 0.0},
    {0, TMatrixDSym(kNPars), allCandidates[1], 1, 0, 0, 0.0, 0.0},
    {0, TMatrixDSym(kNPars), allCandidates[2], 1, 0, 0, 0.0, 0.0},
    {0, TMatrixDSym(kNPars), allCandidates[3], 1, 0, 0, 0.0, 0.0},
    {0, TMatrixDSym(kNPars), allCandidates[4], 1, 0, 0, 0.0, 0.0},
    {0, TMatrixDSym(kNPars), allCandidates[5], 1, 0, 0, 0.0, 0.0},
    {0, TMatrixDSym(kNPars), allCandidates[6], 1, 0, 0, 0.0, 0.0},
    {0, TMatrixDSym(kNPars), allCandidates[7], 1, 0, 0, 0.0, 0.0}
  };
  PhotoPeakFitTrial *bestTrial = 0;
  PhotoPeakFitTrial *fallbackTrial = 0;

  for (int i = 0; i < 8; ++i) {
    PhotoPeakFitTrial &trial = trials[i];
    trial.freePars = PhotoPeakConfiguredFreeParameterCount(trial.candidate, config);
    const bool isBasic = !trial.candidate.useTail &&
      !trial.candidate.useStep && !trial.candidate.useQuadBg;
    const bool selectedByMode =
      mode == kPhotoPeakAuto ||
      (mode == kPhotoPeakHighStat && i == 7) ||
      (mode == kPhotoPeakLowStat && i == 0);
    if (!selectedByMode) {
      continue;
    }
    if (mode == kPhotoPeakAuto && !isBasic && nFitBins <= trial.freePars) {
      continue;
    }

    TString funcName;
    funcName.Form("PhotoPeak_total_fit_%s", trial.candidate.name);
    trial.func = new TF1(funcName.Data(), PhotoPeakEval, fitLow, fitHigh, kNPars);
    PhotoPeakConfigureFunction(trial.func, trial.candidate, a0, b0, c0, r0,
                               beta0, step0, peak0, w0, h0, fitLow, fitHigh,
                               range, hUpper, config);
    TFitResultPtr result = hist->Fit(trial.func, "RQSN");
    trial.status = int(result);
    trial.chi2 = trial.func->GetChisquare();
    trial.ndf = trial.func->GetNDF();
    trial.reducedChi2 = trial.ndf > 0 ? trial.chi2 / trial.ndf :
      std::numeric_limits<double>::quiet_NaN();
    if (result.Get()) {
      trial.cov = result->GetCovarianceMatrix();
    }
    if (isBasic) {
      fallbackTrial = &trial;
    }
    if (mode == kPhotoPeakAuto) {
      if (PhotoPeakTrialIsBetter(&trial, bestTrial)) {
        bestTrial = &trial;
      }
    } else {
      bestTrial = &trial;
    }
  }

  if (!bestTrial) {
    bestTrial = fallbackTrial;
  }
  if (!bestTrial || !bestTrial->func) {
    std::printf("photopeakfit ERROR: no fitting function was available.\n");
    return 3;
  }

  for (int i = 0; i < 8; ++i) {
    if (&trials[i] != bestTrial && trials[i].func) {
      delete trials[i].func;
      trials[i].func = 0;
    }
  }

  TF1 *total = bestTrial->func;
  PhotoPeakWarnInactiveConfig(bestTrial->candidate, config);
  total->SetLineColor(kRed);
  total->SetLineStyle(1);
  total->SetLineWidth(2);
  const int status = bestTrial->status;
  const double chi2 = bestTrial->chi2;
  const int ndf = bestTrial->ndf;
  const double reducedChi2 = bestTrial->reducedChi2;

  double par[kNPars];
  for (int i = 0; i < kNPars; ++i) {
    par[i] = total->GetParameter(i);
  }

  TMatrixDSym cov = bestTrial->cov;
  const double binWidth = hist->GetXaxis()->GetBinWidth(hist->GetXaxis()->FindFixBin(par[kP]));
  const double area = PhotoPeakArea(par, binWidth);
  const double areaErr = PhotoPeakAreaUncertainty(par, cov, binWidth);
  const int printOrder[kNPars] = {kP, kH, kW, kR, kBeta, kStep, kA, kB, kC};

  std::printf("\nphotopeakfit result for %s\n", hist->GetName());
  std::printf("Fit range: [%g, %g], initial peak position: %g\n", fitLow, fitHigh, peak0);
  std::printf("Requested mode: %s\n", PhotoPeakModeName(mode));
  std::printf("Fitting function: %s\n", bestTrial->candidate.name);
  std::printf("Fit bins = %d, free parameters = %d\n", nFitBins, bestTrial->freePars);
  std::printf("Fit status: %d\n", status);
  if (ndf > 0) {
    std::printf("chi2 = %.10g, ndf = %d, reduced chisq = %.10g\n", chi2, ndf, reducedChi2);
  } else {
    std::printf("chi2 = %.10g, ndf = %d, reduced chisq = n/a\n", chi2, ndf);
  }
  std::printf("Photopeak area = %.10g +/- %.10g\n", area, areaErr);
  std::printf("\nParameters:\n");
  for (int i = 0; i < kNPars; ++i) {
    const int ipar = printOrder[i];
    std::printf("  %-5s = % .10g +/- %.10g\n",
                total->GetParName(ipar), total->GetParameter(ipar), total->GetParError(ipar));
  }
  std::printf("\n");

  TF1 *bg = new TF1("PhotoPeak_background_fit", PhotoPeakBgEval, fitLow, fitHigh, kNPars);
  bg->SetParameters(par);
  bg->SetLineColor(kBlack);
  bg->SetLineStyle(2);
  bg->SetLineWidth(2);
  bg->SetNpx(2000);

  if (!gPad) {
    new TCanvas("PhotoPeak_canvas", "PhotoPeak_canvas", 900, 650);
    hist->Draw();
  } else {
    gPad->cd();
    if (!hist->TestBit(TH1::kIsZoomed)) {
      hist->Draw();
    }
  }
  total->Draw("same");
  bg->Draw("same");
  if (gPad) {
    gPad->Modified();
    gPad->Update();
  }

  return status;
}

// ============== photopeakfit ==============
// Purpose: Fit one photopeak with a C++ configuration object.
// Inputs: Histogram pointer, lower x range, upper x range, initial peak x, and config.
// Outputs: Fit status; prints fit results and draws total/background functions.
int photopeakfit(TH1 *hist, double fitLow, double fitHigh, double peak0,
                 const PhotoPeakFitConfig &config)
{
  return PhotoPeakFitRun(hist, fitLow, fitHigh, peak0, config);
}

// ============== photopeakfit ==============
// Purpose: Fit one photopeak with a selected public mode.
// Inputs: Histogram pointer, lower x range, upper x range, initial peak x, and mode.
// Outputs: Fit status; prints fit results and draws total/background functions.
int photopeakfit(TH1 *hist, double fitLow, double fitHigh, double peak0,
                 PhotoPeakFitMode mode = kPhotoPeakAuto)
{
  PhotoPeakFitConfig config;
  config.mode = mode;

  return PhotoPeakFitRun(hist, fitLow, fitHigh, peak0, config);
}

// ============== photopeakfit ==============
// Purpose: Fit one photopeak with a selected public mode integer.
// Inputs: Histogram pointer, lower x range, upper x range, initial peak x, and mode.
// Outputs: Fit status; prints fit results and draws total/background functions.
int photopeakfit(TH1 *hist, double fitLow, double fitHigh, double peak0,
                 int mode)
{
  PhotoPeakFitConfig config;
  config.mode = mode;

  return PhotoPeakFitRun(hist, fitLow, fitHigh, peak0, config);
}

// ============== photopeakfit ==============
// Purpose: Fit one photopeak using a text configuration file.
// Inputs: Histogram pointer, lower x range, upper x range, initial peak x, and config file.
// Outputs: Fit status; prints fit results and draws total/background functions.
int photopeakfit(TH1 *hist, double fitLow, double fitHigh, double peak0,
                 const char *configFile)
{
  PhotoPeakFitConfig config;
  const int readStatus = PhotoPeakReadConfig(configFile, config);
  if (readStatus != 0) {
    return 10 + readStatus;
  }

  return PhotoPeakFitRun(hist, fitLow, fitHigh, peak0, config);
}
