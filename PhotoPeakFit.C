// ROOT macro for RadWare/GF3-style single photopeak fitting.

#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TH1.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TObject.h"
#include "TROOT.h"
#include "TSpectrum.h"
#include "TString.h"
#include "TVirtualPad.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include <initializer_list>

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

enum PhotoPeakTspectrumBgMode {
  kPhotoPeakTspectrumBgNone = 0,
  kPhotoPeakTspectrumBgGlobal = 1,
  kPhotoPeakTspectrumBgLocal = 2
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

struct PhotoPeakIndexedParConfig {
  bool hasInit;
  double init;
  bool hasLimit;
  double limitLow;
  double limitHigh;
  bool hasFix;
  double fix;

  PhotoPeakIndexedParConfig()
  {
    hasInit = false;
    init = 0.0;
    hasLimit = false;
    limitLow = 0.0;
    limitHigh = 0.0;
    hasFix = false;
    fix = 0.0;
  }
};

struct PhotoPeakFitConfig {
  int mode;
  std::string rootFitOption;
  int tspectrumBgMode;
  bool hasTspectrumIteration;
  int tspectrumIteration;
  bool hasTspectrumRange;
  double tspectrumRangeLow;
  double tspectrumRangeHigh;
  bool hasInit[kNPars];
  double init[kNPars];
  bool hasLimit[kNPars];
  double limitLow[kNPars];
  double limitHigh[kNPars];
  bool hasFix[kNPars];
  double fix[kNPars];
  bool relativePosition;
  bool relativeFwhm;
  bool hasWScaleInit;
  double wScaleInit;
  bool hasWScaleLimit;
  double wScaleLimitLow;
  double wScaleLimitHigh;
  bool hasWScaleFix;
  double wScaleFix;
  std::vector<PhotoPeakIndexedParConfig> multiP;
  std::vector<PhotoPeakIndexedParConfig> multiW;
  std::vector<PhotoPeakIndexedParConfig> multiH;

  PhotoPeakFitConfig()
  {
    mode = kPhotoPeakAuto;
    rootFitOption = "RQSN";
    tspectrumBgMode = kPhotoPeakTspectrumBgNone;
    hasTspectrumIteration = false;
    tspectrumIteration = 20;
    hasTspectrumRange = false;
    tspectrumRangeLow = 0.0;
    tspectrumRangeHigh = 0.0;
    for (int i = 0; i < kNPars; ++i) {
      hasInit[i] = false;
      init[i] = 0.0;
      hasLimit[i] = false;
      limitLow[i] = 0.0;
      limitHigh[i] = 0.0;
      hasFix[i] = false;
      fix[i] = 0.0;
    }
    relativePosition = false;
    relativeFwhm = true;
    hasWScaleInit = false;
    wScaleInit = 1.0;
    hasWScaleLimit = false;
    wScaleLimitLow = 0.0;
    wScaleLimitHigh = 0.0;
    hasWScaleFix = false;
    wScaleFix = 1.0;
  }
};

struct PhotoPeakMultiParMap {
  int nPeaks;
  bool relativePosition;
  bool relativeFwhm;
  bool useTail;
  bool useStep;
  bool useQuadBg;
  bool hasFixedPositionAnchor;
  int fixedPositionAnchor;
  int wScaleIndex;
  std::vector<int> pIndex;
  std::vector<int> wIndex;
  std::vector<int> hIndex;
  std::vector<double> p0;
  std::vector<double> w0;
  std::vector<double> pFixed;
  std::vector<double> wFixed;
  std::vector<bool> hasPFixed;
  std::vector<bool> hasWFixed;
  std::vector<bool> hasHFixed;
  std::vector<std::string> parNames;

  PhotoPeakMultiParMap()
  {
    nPeaks = 0;
    relativePosition = false;
    relativeFwhm = true;
    useTail = false;
    useStep = false;
    useQuadBg = false;
    hasFixedPositionAnchor = false;
    fixedPositionAnchor = -1;
    wScaleIndex = -1;
  }
};

struct PhotoPeakDerivedValue {
  double value;
  double error;
  const char *source;

  PhotoPeakDerivedValue()
  {
    value = 0.0;
    error = 0.0;
    source = "fit";
  }
};

double gPhotoPeakFitLow = 0.0;
double gPhotoPeakFitHigh = 0.0;
TH1 *gPhotoPeakTspectrumBg = 0;
PhotoPeakMultiParMap *gPhotoPeakMultiMap = 0;
PhotoPeakMultiParMap gPhotoPeakMultiMapStorage;

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

// ============== PhotoPeakStripQuotes ==============
// Purpose: Remove one matching pair of surrounding quotes from a string.
// Inputs: String value.
// Outputs: Unquoted string when quotes are present.
std::string PhotoPeakStripQuotes(const std::string &text)
{
  const std::string trimmed = PhotoPeakTrim(text);
  if (trimmed.size() >= 2 &&
      ((trimmed[0] == '"' && trimmed[trimmed.size() - 1] == '"') ||
       (trimmed[0] == '\'' && trimmed[trimmed.size() - 1] == '\''))) {
    return trimmed.substr(1, trimmed.size() - 2);
  }

  return trimmed;
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

// ============== PhotoPeakBoolName ==============
// Purpose: Convert a boolean value to a printable true/false string.
// Inputs: Boolean value.
// Outputs: Static boolean string.
const char *PhotoPeakBoolName(bool value)
{
  return value ? "true" : "false";
}

// ============== PhotoPeakParseBool ==============
// Purpose: Convert a user string to a boolean value.
// Inputs: Text and output boolean reference.
// Outputs: True when the text is recognized.
bool PhotoPeakParseBool(const std::string &text, bool &value)
{
  const std::string key = PhotoPeakLower(PhotoPeakTrim(text));
  if (key == "true" || key == "yes" || key == "on" || key == "1") {
    value = true;
    return true;
  }
  if (key == "false" || key == "no" || key == "off" || key == "0") {
    value = false;
    return true;
  }

  return false;
}

// ============== PhotoPeakParseIndexedParName ==============
// Purpose: Parse P[i], W[i], H[i], or WSCALE parameter names.
// Inputs: Parameter name plus output scalar/index references.
// Outputs: True when the parameter name is a multi-peak parameter.
bool PhotoPeakParseIndexedParName(const std::string &name,
                                  int &ipar, int &ipeak, bool &isWScale)
{
  const std::string key = PhotoPeakUpper(PhotoPeakTrim(name));
  ipar = -1;
  ipeak = -1;
  isWScale = false;
  if (key == "WSCALE" || key == "W_SCALE" || key == "FWHMSCALE" ||
      key == "FWHM_SCALE") {
    isWScale = true;
    return true;
  }
  if (key.size() < 4 || key[1] != '[' || key[key.size() - 1] != ']') {
    return false;
  }
  if (key[0] == 'P') {
    ipar = kP;
  } else if (key[0] == 'W') {
    ipar = kW;
  } else if (key[0] == 'H') {
    ipar = kH;
  } else {
    return false;
  }

  const std::string indexText = key.substr(2, key.size() - 3);
  if (indexText.empty()) {
    return false;
  }
  for (size_t i = 0; i < indexText.size(); ++i) {
    if (!std::isdigit(static_cast<unsigned char>(indexText[i]))) {
      return false;
    }
  }
  ipeak = std::atoi(indexText.c_str());

  return ipeak >= 0;
}

// ============== PhotoPeakEnsureIndexedPar ==============
// Purpose: Ensure indexed multi-peak config storage exists.
// Inputs: Config vector and peak index.
// Outputs: Reference to the requested indexed setting.
PhotoPeakIndexedParConfig &PhotoPeakEnsureIndexedPar(std::vector<PhotoPeakIndexedParConfig> &settings,
                                                     int ipeak)
{
  if (static_cast<int>(settings.size()) <= ipeak) {
    settings.resize(ipeak + 1);
  }

  return settings[ipeak];
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

// ============== PhotoPeakTspectrumBgValue ==============
// Purpose: Read the TSpectrum background estimate at one x-value.
// Inputs: X-value.
// Outputs: Background bin content, or zero when no estimate exists.
double PhotoPeakTspectrumBgValue(double x)
{
  if (!gPhotoPeakTspectrumBg) {
    return 0.0;
  }
  const int bin = gPhotoPeakTspectrumBg->GetXaxis()->FindFixBin(x);
  const int clampedBin = std::max(1, std::min(gPhotoPeakTspectrumBg->GetNbinsX(), bin));

  return gPhotoPeakTspectrumBg->GetBinContent(clampedBin);
}

// ============== PhotoPeakDisplayEval ==============
// Purpose: Evaluate the total fit on the original histogram scale.
// Inputs: ROOT x array and parameter array.
// Outputs: Fit value plus removed TSpectrum background.
double PhotoPeakDisplayEval(double *x, double *par)
{
  return PhotoPeakEval(x, par) + PhotoPeakTspectrumBgValue(x[0]);
}

// ============== PhotoPeakBgDisplayEval ==============
// Purpose: Evaluate the background fit on the original histogram scale.
// Inputs: ROOT x array and parameter array.
// Outputs: Residual fit background plus removed TSpectrum background.
double PhotoPeakBgDisplayEval(double *x, double *par)
{
  return PhotoPeakBgEval(x, par) + PhotoPeakTspectrumBgValue(x[0]);
}

// ============== PhotoPeakTspectrumBgEval ==============
// Purpose: Evaluate only the removed TSpectrum background estimate.
// Inputs: ROOT x array and unused parameter array.
// Outputs: Background estimate at x.
double PhotoPeakTspectrumBgEval(double *x, double *)
{
  return PhotoPeakTspectrumBgValue(x[0]);
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

struct PhotoPeakTspectrumBgResult {
  TH1 *fitHist;
  TH1 *bgHist;
  bool ownsFitHist;
  int mode;
  int iteration;
  double rangeLow;
  double rangeHigh;

  PhotoPeakTspectrumBgResult()
  {
    fitHist = 0;
    bgHist = 0;
    ownsFitHist = false;
    mode = kPhotoPeakTspectrumBgNone;
    iteration = 20;
    rangeLow = 0.0;
    rangeHigh = 0.0;
  }
};

// ============== PhotoPeakClampRangeToHist ==============
// Purpose: Clamp an x range to the histogram axis limits.
// Inputs: Histogram pointer and x range references.
// Outputs: Clamped x range.
void PhotoPeakClampRangeToHist(TH1 *hist, double &low, double &high)
{
  const double xMin = hist->GetXaxis()->GetXmin();
  const double xMax = hist->GetXaxis()->GetXmax();
  low = std::max(xMin, std::min(xMax, low));
  high = std::max(xMin, std::min(xMax, high));
  if (low > high) {
    std::swap(low, high);
  }
}

// ============== PhotoPeakReferenceFwhm ==============
// Purpose: Choose the FWHM used for TSpectrum local defaults.
// Inputs: Initial peak x-value and fit configuration.
// Outputs: Positive FWHM estimate.
double PhotoPeakReferenceFwhm(double peak0, const PhotoPeakFitConfig &config)
{
  if (config.hasFix[kW] && config.fix[kW] > 0.0) {
    return config.fix[kW];
  }

  return std::sqrt(std::max(9.0 + 0.004 * peak0, 1.0e-12));
}

// ============== PhotoPeakDefaultLocalIteration ==============
// Purpose: Calculate the automatic local TSpectrum iteration count.
// Inputs: FWHM estimate and histogram bin width.
// Outputs: Clamped iteration count.
int PhotoPeakDefaultLocalIteration(double fwhm, double binWidth)
{
  const double fwhmBins = fwhm / std::max(binWidth, 1.0e-12);
  const int estimate = static_cast<int>(std::lround(2.5 * fwhmBins));

  return std::max(6, std::min(20, estimate));
}

// ============== PhotoPeakPrepareTspectrumBackground ==============
// Purpose: Build a background-subtracted fit histogram when requested.
// Inputs: Histogram, fit limits, peak estimate, and configuration.
// Outputs: Fit histogram plus TSpectrum background metadata.
PhotoPeakTspectrumBgResult PhotoPeakPrepareTspectrumBackground(TH1 *hist,
                                                               double fitLow,
                                                               double fitHigh,
                                                               double peak0,
                                                               const PhotoPeakFitConfig &config)
{
  PhotoPeakTspectrumBgResult result;
  result.fitHist = hist;
  result.mode = config.tspectrumBgMode;
  if (result.mode != kPhotoPeakTspectrumBgNone &&
      result.mode != kPhotoPeakTspectrumBgGlobal &&
      result.mode != kPhotoPeakTspectrumBgLocal) {
    result.mode = kPhotoPeakTspectrumBgNone;
  }
  if (result.mode == kPhotoPeakTspectrumBgNone) {
    return result;
  }

  const int peakBin = hist->GetXaxis()->FindFixBin(peak0);
  const double binWidth = hist->GetXaxis()->GetBinWidth(peakBin);
  const double fwhm = PhotoPeakReferenceFwhm(peak0, config);
  if (result.mode == kPhotoPeakTspectrumBgGlobal) {
    result.iteration = config.hasTspectrumIteration ? config.tspectrumIteration : 20;
    result.rangeLow = hist->GetXaxis()->GetXmin();
    result.rangeHigh = hist->GetXaxis()->GetXmax();
  } else {
    result.iteration = config.hasTspectrumIteration ?
      config.tspectrumIteration : PhotoPeakDefaultLocalIteration(fwhm, binWidth);
    if (config.hasTspectrumRange) {
      result.rangeLow = config.tspectrumRangeLow;
      result.rangeHigh = config.tspectrumRangeHigh;
    } else {
      const double fitWidth = fitHigh - fitLow;
      const double padding = std::max(std::max(5.0 * fwhm, 0.5 * fitWidth),
                                      result.iteration * binWidth);
      result.rangeLow = fitLow - padding;
      result.rangeHigh = fitHigh + padding;
    }
    if (result.rangeLow > fitLow || result.rangeHigh < fitHigh) {
      std::printf("photopeakfit WARNING: local TSpectrum range does not cover fit range; expanding it.\n");
      result.rangeLow = std::min(result.rangeLow, fitLow);
      result.rangeHigh = std::max(result.rangeHigh, fitHigh);
    }
    PhotoPeakClampRangeToHist(hist, result.rangeLow, result.rangeHigh);
  }

  TString sourceName;
  sourceName.Form("%s_tspectrum_source", hist->GetName());
  TH1 *source = static_cast<TH1 *>(hist->Clone(sourceName.Data()));
  source->SetDirectory(0);
  if (result.mode == kPhotoPeakTspectrumBgGlobal) {
    source->GetXaxis()->SetRange(1, source->GetNbinsX());
  } else if (result.mode == kPhotoPeakTspectrumBgLocal) {
    const int lowBin = source->GetXaxis()->FindFixBin(result.rangeLow);
    const int highBin = source->GetXaxis()->FindFixBin(result.rangeHigh);
    source->GetXaxis()->SetRange(lowBin, highBin);
  }

  TSpectrum spectrum;
  result.bgHist = spectrum.Background(source, result.iteration, "");
  delete source;
  if (!result.bgHist) {
    std::printf("photopeakfit WARNING: TSpectrum background failed; fitting original histogram.\n");
    result.mode = kPhotoPeakTspectrumBgNone;
    result.fitHist = hist;
    return result;
  }
  result.bgHist->SetDirectory(0);
  result.bgHist->SetName("PhotoPeak_tspectrum_background");

  TString fitName;
  fitName.Form("%s_tspectrum_subtracted", hist->GetName());
  result.fitHist = static_cast<TH1 *>(hist->Clone(fitName.Data()));
  result.fitHist->SetDirectory(0);
  result.ownsFitHist = true;

  int subtractLowBin = 1;
  int subtractHighBin = hist->GetNbinsX();
  if (result.mode == kPhotoPeakTspectrumBgLocal) {
    subtractLowBin = hist->GetXaxis()->FindFixBin(fitLow);
    subtractHighBin = hist->GetXaxis()->FindFixBin(fitHigh);
    if (subtractLowBin > subtractHighBin) {
      std::swap(subtractLowBin, subtractHighBin);
    }
  }
  subtractLowBin = std::max(1, std::min(hist->GetNbinsX(), subtractLowBin));
  subtractHighBin = std::max(1, std::min(hist->GetNbinsX(), subtractHighBin));

  for (int bin = subtractLowBin; bin <= subtractHighBin; ++bin) {
    const double value = hist->GetBinContent(bin) - result.bgHist->GetBinContent(bin);
    result.fitHist->SetBinContent(bin, value);
  }

  return result;
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

// ============== PhotoPeakTspectrumBgModeName ==============
// Purpose: Convert a TSpectrum background mode value to a printable name.
// Inputs: Background mode integer.
// Outputs: Static mode name string.
const char *PhotoPeakTspectrumBgModeName(int mode)
{
  if (mode == kPhotoPeakTspectrumBgGlobal) {
    return "global";
  }
  if (mode == kPhotoPeakTspectrumBgLocal) {
    return "local";
  }

  return "none";
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

// ============== PhotoPeakParseTspectrumBgMode ==============
// Purpose: Convert a user string to a TSpectrum background mode.
// Inputs: Mode string and output mode reference.
// Outputs: True when the mode string is recognized.
bool PhotoPeakParseTspectrumBgMode(const std::string &text, int &mode)
{
  const std::string key = PhotoPeakLower(PhotoPeakTrim(text));
  if (key == "none" || key == "off" || key == "no") {
    mode = kPhotoPeakTspectrumBgNone;
    return true;
  }
  if (key == "global") {
    mode = kPhotoPeakTspectrumBgGlobal;
    return true;
  }
  if (key == "local") {
    mode = kPhotoPeakTspectrumBgLocal;
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
      if (line[i] == '=' || line[i] == ',') {
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

    if (keyword == "rootfitoption" || keyword == "rootoption" ||
        keyword == "fitoption") {
      std::string optionText;
      words >> optionText;
      optionText = PhotoPeakStripQuotes(optionText);
      if (optionText.empty()) {
        std::printf("photopeakfit WARNING: %s:%d empty ROOT fit option; using RQSN.\n",
                    configFile, lineNumber);
        optionText = "RQSN";
      }
      config.rootFitOption = optionText;
      continue;
    }

    if (keyword == "option") {
      std::string optionText;
      words >> optionText;
      int parsedMode = kPhotoPeakTspectrumBgNone;
      if (!PhotoPeakParseTspectrumBgMode(optionText, parsedMode)) {
        std::printf("photopeakfit WARNING: %s:%d unknown TSpectrum background option '%s'; using none.\n",
                    configFile, lineNumber, optionText.c_str());
        parsedMode = kPhotoPeakTspectrumBgNone;
      }
      config.tspectrumBgMode = parsedMode;
      continue;
    }

    if (keyword == "iteration" || keyword == "iterations") {
      int iteration = 20;
      if (!(words >> iteration)) {
        std::printf("photopeakfit WARNING: %s:%d missing TSpectrum iteration; using default.\n",
                    configFile, lineNumber);
        continue;
      }
      config.hasTspectrumIteration = true;
      config.tspectrumIteration = std::max(1, iteration);
      continue;
    }

    if (keyword == "range") {
      double low = 0.0;
      double high = 0.0;
      if (!(words >> low >> high)) {
        std::printf("photopeakfit WARNING: %s:%d missing TSpectrum range values; ignoring line.\n",
                    configFile, lineNumber);
        continue;
      }
      if (low > high) {
        std::swap(low, high);
      }
      config.hasTspectrumRange = true;
      config.tspectrumRangeLow = low;
      config.tspectrumRangeHigh = high;
      continue;
    }

    if (keyword == "relativeposition" || keyword == "relativepos" ||
        keyword == "relpos") {
      std::string valueText;
      words >> valueText;
      bool value = false;
      if (!PhotoPeakParseBool(valueText, value)) {
        std::printf("photopeakfit WARNING: %s:%d unknown relativePosition value '%s'; using false.\n",
                    configFile, lineNumber, valueText.c_str());
        value = false;
      }
      config.relativePosition = value;
      continue;
    }

    if (keyword == "relativefwhm" || keyword == "relativewidth" ||
        keyword == "relwid" || keyword == "relfwhm") {
      std::string valueText;
      words >> valueText;
      bool value = true;
      if (!PhotoPeakParseBool(valueText, value)) {
        std::printf("photopeakfit WARNING: %s:%d unknown relativeFwhm value '%s'; using true.\n",
                    configFile, lineNumber, valueText.c_str());
        value = true;
      }
      config.relativeFwhm = value;
      continue;
    }

    if (keyword == "init" || keyword == "fix" || keyword == "limit") {
      std::string parName;
      words >> parName;
      int indexedPar = -1;
      int indexedPeak = -1;
      bool isWScale = false;
      if (PhotoPeakParseIndexedParName(parName, indexedPar, indexedPeak, isWScale)) {
        if (keyword == "init" || keyword == "fix") {
          double value = 0.0;
          if (!(words >> value)) {
            std::printf("photopeakfit WARNING: %s:%d missing value; ignoring line.\n",
                        configFile, lineNumber);
            continue;
          }
          if (isWScale) {
            if (keyword == "init") {
              config.hasWScaleInit = true;
              config.wScaleInit = value;
            } else {
              config.hasWScaleFix = true;
              config.wScaleFix = value;
            }
            continue;
          }
          std::vector<PhotoPeakIndexedParConfig> *settings = 0;
          if (indexedPar == kP) {
            settings = &config.multiP;
          } else if (indexedPar == kW) {
            settings = &config.multiW;
          } else if (indexedPar == kH) {
            settings = &config.multiH;
          }
          PhotoPeakIndexedParConfig &entry =
            PhotoPeakEnsureIndexedPar(*settings, indexedPeak);
          if (keyword == "init") {
            entry.hasInit = true;
            entry.init = value;
          } else {
            entry.hasFix = true;
            entry.fix = value;
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
        if (isWScale) {
          config.hasWScaleLimit = true;
          config.wScaleLimitLow = low;
          config.wScaleLimitHigh = high;
          continue;
        }
        std::vector<PhotoPeakIndexedParConfig> *settings = 0;
        if (indexedPar == kP) {
          settings = &config.multiP;
        } else if (indexedPar == kW) {
          settings = &config.multiW;
        } else if (indexedPar == kH) {
          settings = &config.multiH;
        }
        PhotoPeakIndexedParConfig &entry =
          PhotoPeakEnsureIndexedPar(*settings, indexedPeak);
        entry.hasLimit = true;
        entry.limitLow = low;
        entry.limitHigh = high;
        continue;
      }

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
    if (config.hasInit[ipar] && ipar != kP) {
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

// ============== PhotoPeakIndexedEntry ==============
// Purpose: Read an optional indexed multi-peak configuration entry.
// Inputs: Config vector and peak index.
// Outputs: Pointer to the entry, or null when absent.
const PhotoPeakIndexedParConfig *PhotoPeakIndexedEntry(const std::vector<PhotoPeakIndexedParConfig> &settings,
                                                       int ipeak)
{
  if (ipeak < 0 || ipeak >= static_cast<int>(settings.size())) {
    return 0;
  }

  return &settings[ipeak];
}

// ============== PhotoPeakPeakShape ==============
// Purpose: Evaluate one GF3-style photopeak component.
// Inputs: X-value, shared shape parameters, peak position, FWHM, and height.
// Outputs: Photopeak value without background.
double PhotoPeakPeakShape(double xx, double rPercent, double betaValue,
                          double p, double wFwhm, double h)
{
  const double beta = std::max(betaValue, 1.0e-12);
  const double fwhm = std::max(wFwhm, 1.0e-12);
  const double sigma = fwhm / 2.35482;
  const double w = (xx - p) / (sigma * TMath::Sqrt2());
  const double y = fwhm / (beta * 3.33021838);
  const double erfcY = std::max(TMath::Erfc(y), 1.0e-300);
  const double gaussian = std::exp(-w * w);
  double skew = 0.0;
  const double tailArg = (xx - p) / beta;
  if (std::fabs(tailArg) <= 700.0) {
    skew = std::exp(tailArg) * TMath::Erfc(w + y) / erfcY;
  }
  const double r = rPercent / 100.0;

  return h * ((1.0 - r) * gaussian + r * skew);
}

// ============== PhotoPeakMultiPositionValue ==============
// Purpose: Calculate one multi-peak centroid from the active parameter map.
// Inputs: Peak index and fit parameter array.
// Outputs: Centroid value.
double PhotoPeakMultiPositionValue(int ipeak, const double *par)
{
  const PhotoPeakMultiParMap &map = *gPhotoPeakMultiMap;
  if (map.hasPFixed[ipeak]) {
    return map.pFixed[ipeak];
  }
  if (!map.relativePosition) {
    return par[map.pIndex[ipeak]];
  }
  if (map.hasFixedPositionAnchor) {
    const int anchor = map.fixedPositionAnchor;
    return map.pFixed[anchor] + (map.p0[ipeak] - map.p0[anchor]);
  }

  return par[map.pIndex[0]] + (map.p0[ipeak] - map.p0[0]);
}

// ============== PhotoPeakMultiFwhmValue ==============
// Purpose: Calculate one multi-peak FWHM from the active parameter map.
// Inputs: Peak index and fit parameter array.
// Outputs: FWHM value.
double PhotoPeakMultiFwhmValue(int ipeak, const double *par)
{
  const PhotoPeakMultiParMap &map = *gPhotoPeakMultiMap;
  if (map.hasWFixed[ipeak]) {
    return map.wFixed[ipeak];
  }
  if (!map.relativeFwhm) {
    return par[map.wIndex[ipeak]];
  }
  if (map.wScaleIndex >= 0) {
    return par[map.wScaleIndex] * map.w0[ipeak];
  }

  return map.w0[ipeak];
}

// ============== PhotoPeakMultiPeakHeight ==============
// Purpose: Read one multi-peak height from the active parameter map.
// Inputs: Peak index and fit parameter array.
// Outputs: Height value.
double PhotoPeakMultiPeakHeight(int ipeak, const double *par)
{
  const PhotoPeakMultiParMap &map = *gPhotoPeakMultiMap;
  return par[map.hIndex[ipeak]];
}

// ============== PhotoPeakMultiEval ==============
// Purpose: Evaluate the total multi-peak GF3-style fitting function.
// Inputs: ROOT x array and parameter array.
// Outputs: Total multi-peak fit value at x.
double PhotoPeakMultiEval(double *x, double *par)
{
  const PhotoPeakMultiParMap &map = *gPhotoPeakMultiMap;
  const double xx = x[0];
  const double xCentered = xx - 0.5 * (gPhotoPeakFitLow + gPhotoPeakFitHigh);
  double value = par[kA] + par[kB] * xCentered + par[kC] * xCentered * xCentered;
  for (int ipeak = 0; ipeak < map.nPeaks; ++ipeak) {
    const double p = PhotoPeakMultiPositionValue(ipeak, par);
    const double wFwhm = PhotoPeakMultiFwhmValue(ipeak, par);
    const double h = par[map.hIndex[ipeak]];
    value += PhotoPeakPeakShape(xx, par[kR], par[kBeta], p, wFwhm, h);
    const double sigma = std::max(wFwhm, 1.0e-12) / 2.35482;
    const double w = (xx - p) / (sigma * TMath::Sqrt2());
    value += h * par[kStep] * TMath::Erfc(w) / 200.0;
  }

  return value;
}

// ============== PhotoPeakMultiBgEval ==============
// Purpose: Evaluate the multi-peak background-only part.
// Inputs: ROOT x array and parameter array.
// Outputs: Shared background plus per-peak step background.
double PhotoPeakMultiBgEval(double *x, double *par)
{
  const PhotoPeakMultiParMap &map = *gPhotoPeakMultiMap;
  const double xx = x[0];
  const double xCentered = xx - 0.5 * (gPhotoPeakFitLow + gPhotoPeakFitHigh);
  double value = par[kA] + par[kB] * xCentered + par[kC] * xCentered * xCentered;
  for (int ipeak = 0; ipeak < map.nPeaks; ++ipeak) {
    const double p = PhotoPeakMultiPositionValue(ipeak, par);
    const double wFwhm = PhotoPeakMultiFwhmValue(ipeak, par);
    const double h = par[map.hIndex[ipeak]];
    const double sigma = std::max(wFwhm, 1.0e-12) / 2.35482;
    const double w = (xx - p) / (sigma * TMath::Sqrt2());
    value += h * par[kStep] * TMath::Erfc(w) / 200.0;
  }

  return value;
}

// ============== PhotoPeakMultiDisplayEval ==============
// Purpose: Evaluate the multi-peak total fit on the original histogram scale.
// Inputs: ROOT x array and parameter array.
// Outputs: Fit value plus removed TSpectrum background.
double PhotoPeakMultiDisplayEval(double *x, double *par)
{
  return PhotoPeakMultiEval(x, par) + PhotoPeakTspectrumBgValue(x[0]);
}

// ============== PhotoPeakMultiBgDisplayEval ==============
// Purpose: Evaluate the multi-peak background on the original histogram scale.
// Inputs: ROOT x array and parameter array.
// Outputs: Background value plus removed TSpectrum background.
double PhotoPeakMultiBgDisplayEval(double *x, double *par)
{
  return PhotoPeakMultiBgEval(x, par) + PhotoPeakTspectrumBgValue(x[0]);
}

// ============== PhotoPeakMultiPeakEval ==============
// Purpose: Evaluate one single peak component for drawing.
// Inputs: ROOT x array and parameter array.
// Outputs: One peak component on the original histogram scale baseline.
double PhotoPeakMultiPeakEval(double *x, double *par)
{
  if (!gPhotoPeakMultiMap) {
    return 0.0;
  }
  const int ipeak =
    static_cast<int>(std::lround(par[gPhotoPeakMultiMap->parNames.size()]));
  if (ipeak < 0 || ipeak >= gPhotoPeakMultiMap->nPeaks) {
    return 0.0;
  }
  const double p = PhotoPeakMultiPositionValue(ipeak, par);
  const double wFwhm = PhotoPeakMultiFwhmValue(ipeak, par);
  const double h = par[gPhotoPeakMultiMap->hIndex[ipeak]];

  return PhotoPeakPeakShape(x[0], par[kR], par[kBeta], p, wFwhm, h) +
    PhotoPeakMultiBgEval(x, par) + PhotoPeakTspectrumBgValue(x[0]);
}

// ============== PhotoPeakBuildMultiParMap ==============
// Purpose: Build the logical-to-TF1 parameter map for multi-peak fitting.
// Inputs: Initial peaks, reference FWHM values, candidate, and configuration.
// Outputs: Filled parameter map.
PhotoPeakMultiParMap PhotoPeakBuildMultiParMap(const std::vector<double> &peaks,
                                               const std::vector<double> &w0,
                                               const PhotoPeakFitCandidate &candidate,
                                               const PhotoPeakFitConfig &config)
{
  PhotoPeakMultiParMap map;
  map.nPeaks = static_cast<int>(peaks.size());
  map.relativePosition = config.relativePosition;
  map.relativeFwhm = config.relativeFwhm;
  map.useTail = candidate.useTail;
  map.useStep = candidate.useStep;
  map.useQuadBg = candidate.useQuadBg;
  map.p0 = peaks;
  map.w0 = w0;
  map.pIndex.assign(map.nPeaks, -1);
  map.wIndex.assign(map.nPeaks, -1);
  map.hIndex.assign(map.nPeaks, -1);
  map.pFixed.assign(map.nPeaks, 0.0);
  map.wFixed.assign(map.nPeaks, 0.0);
  map.hasPFixed.assign(map.nPeaks, false);
  map.hasWFixed.assign(map.nPeaks, false);
  map.hasHFixed.assign(map.nPeaks, false);
  map.parNames.push_back("A");
  map.parNames.push_back("B");
  map.parNames.push_back("C");
  map.parNames.push_back("R");
  map.parNames.push_back("BETA");
  map.parNames.push_back("STEP");

  for (int ipeak = 0; ipeak < map.nPeaks; ++ipeak) {
    const PhotoPeakIndexedParConfig *pEntry =
      PhotoPeakIndexedEntry(config.multiP, ipeak);
    const PhotoPeakIndexedParConfig *wEntry =
      PhotoPeakIndexedEntry(config.multiW, ipeak);
    const PhotoPeakIndexedParConfig *hEntry =
      PhotoPeakIndexedEntry(config.multiH, ipeak);
    if (pEntry && pEntry->hasFix) {
      map.hasPFixed[ipeak] = true;
      map.pFixed[ipeak] = pEntry->fix;
      if (!map.hasFixedPositionAnchor) {
        map.hasFixedPositionAnchor = true;
        map.fixedPositionAnchor = ipeak;
      }
    }
    if (wEntry && wEntry->hasFix) {
      map.hasWFixed[ipeak] = true;
      map.wFixed[ipeak] = wEntry->fix;
    }
    if (hEntry && hEntry->hasFix) {
      map.hasHFixed[ipeak] = true;
    }
  }

  if (map.relativePosition) {
    if (!map.hasFixedPositionAnchor) {
      map.pIndex[0] = static_cast<int>(map.parNames.size());
      map.parNames.push_back("PMASTER");
    }
  } else {
    for (int ipeak = 0; ipeak < map.nPeaks; ++ipeak) {
      if (!map.hasPFixed[ipeak]) {
        map.pIndex[ipeak] = static_cast<int>(map.parNames.size());
        TString name;
        name.Form("P[%d]", ipeak);
        map.parNames.push_back(name.Data());
      }
    }
  }

  if (map.relativeFwhm) {
    bool needsWScale = false;
    for (int ipeak = 0; ipeak < map.nPeaks; ++ipeak) {
      if (!map.hasWFixed[ipeak]) {
        needsWScale = true;
      }
    }
    if (needsWScale) {
      map.wScaleIndex = static_cast<int>(map.parNames.size());
      map.parNames.push_back("WSCALE");
    }
  } else {
    for (int ipeak = 0; ipeak < map.nPeaks; ++ipeak) {
      if (!map.hasWFixed[ipeak]) {
        map.wIndex[ipeak] = static_cast<int>(map.parNames.size());
        TString name;
        name.Form("W[%d]", ipeak);
        map.parNames.push_back(name.Data());
      }
    }
  }

  for (int ipeak = 0; ipeak < map.nPeaks; ++ipeak) {
    map.hIndex[ipeak] = static_cast<int>(map.parNames.size());
    TString name;
    name.Form("H[%d]", ipeak);
    map.parNames.push_back(name.Data());
  }

  return map;
}

// ============== PhotoPeakMultiFreeParameterCount ==============
// Purpose: Count free parameters in a configured multi-peak map.
// Inputs: Parameter map, candidate, and configuration.
// Outputs: Number of free parameters.
int PhotoPeakMultiFreeParameterCount(const PhotoPeakMultiParMap &map,
                                     const PhotoPeakFitCandidate &candidate,
                                     const PhotoPeakFitConfig &config)
{
  int freePars = static_cast<int>(map.parNames.size());
  if (!candidate.useQuadBg || config.hasFix[kC]) {
    --freePars;
  }
  if (!candidate.useTail || (config.hasFix[kR] && config.fix[kR] == 0.0)) {
    freePars -= 2;
  } else {
    if (config.hasFix[kR]) {
      --freePars;
    }
    if (config.hasFix[kBeta]) {
      --freePars;
    }
  }
  if (!candidate.useStep || config.hasFix[kStep]) {
    --freePars;
  }
  if (config.hasFix[kA]) {
    --freePars;
  }
  if (config.hasFix[kB]) {
    --freePars;
  }
  if (map.wScaleIndex >= 0 && config.hasWScaleFix) {
    --freePars;
  }
  for (int ipeak = 0; ipeak < map.nPeaks; ++ipeak) {
    if (map.hasHFixed[ipeak]) {
      --freePars;
    }
  }

  return std::max(0, freePars);
}

// ============== PhotoPeakSetMultiParNames ==============
// Purpose: Assign TF1 parameter names from the multi-peak map.
// Inputs: TF1 pointer and parameter map.
// Outputs: Named TF1 parameters.
void PhotoPeakSetMultiParNames(TF1 *func, const PhotoPeakMultiParMap &map)
{
  for (int ipar = 0; ipar < static_cast<int>(map.parNames.size()); ++ipar) {
    func->SetParName(ipar, map.parNames[ipar].c_str());
  }
}

// ============== PhotoPeakConfigureMultiFunction ==============
// Purpose: Set initial values, limits, and fixed parameters for multi-peak fitting.
// Inputs: Function, map, candidate, defaults, range, and config.
// Outputs: Configured ROOT function.
void PhotoPeakConfigureMultiFunction(TF1 *func, const PhotoPeakMultiParMap &map,
                                     const PhotoPeakFitCandidate &candidate,
                                     double a0, double b0, double c0,
                                     double r0, double beta0, double step0,
                                     const std::vector<double> &h0,
                                     double fitLow, double fitHigh,
                                     double range, double hUpper,
                                     const PhotoPeakFitConfig &config)
{
  PhotoPeakSetMultiParNames(func, map);
  func->SetParameter(kA, config.hasInit[kA] ? config.init[kA] : a0);
  func->SetParameter(kB, config.hasInit[kB] ? config.init[kB] : b0);
  func->SetParameter(kC, config.hasInit[kC] ? config.init[kC] : c0);
  func->SetParameter(kR, config.hasInit[kR] ? config.init[kR] : r0);
  func->SetParameter(kBeta, config.hasInit[kBeta] ? config.init[kBeta] : beta0);
  func->SetParameter(kStep, config.hasInit[kStep] ? config.init[kStep] : step0);
  func->SetParLimits(kR, 0.0, 100.0);
  func->SetParLimits(kBeta, 1.0e-6, 10.0 * range);
  func->SetParLimits(kStep, 0.0, 100.0);
  if (!candidate.useQuadBg) {
    func->FixParameter(kC, 0.0);
  }
  if (!candidate.useTail || (config.hasFix[kR] && config.fix[kR] == 0.0)) {
    func->FixParameter(kR, 0.0);
    func->FixParameter(kBeta, config.hasInit[kBeta] ? config.init[kBeta] : beta0);
  }
  if (!candidate.useStep) {
    func->FixParameter(kStep, 0.0);
  }
  for (int ipar = kA; ipar <= kStep; ++ipar) {
    if ((ipar == kC && !candidate.useQuadBg) ||
        ((ipar == kR || ipar == kBeta) && !candidate.useTail) ||
        (ipar == kStep && !candidate.useStep)) {
      continue;
    }
    if (config.hasLimit[ipar]) {
      func->SetParLimits(ipar, config.limitLow[ipar], config.limitHigh[ipar]);
    }
    if (config.hasFix[ipar]) {
      if (!(ipar == kBeta && config.hasFix[kR] && config.fix[kR] == 0.0)) {
        func->FixParameter(ipar, config.fix[ipar]);
      }
    }
  }

  if (map.relativePosition) {
    if (!map.hasFixedPositionAnchor) {
      const int pIndex = map.pIndex[0];
      func->SetParameter(pIndex, map.p0[0]);
      func->SetParLimits(pIndex, fitLow, fitHigh);
      const PhotoPeakIndexedParConfig *entry =
        PhotoPeakIndexedEntry(config.multiP, 0);
      if (entry && entry->hasLimit) {
        func->SetParLimits(pIndex, entry->limitLow, entry->limitHigh);
      }
    }
  } else {
    for (int ipeak = 0; ipeak < map.nPeaks; ++ipeak) {
      if (map.hasPFixed[ipeak]) {
        continue;
      }
      const int pIndex = map.pIndex[ipeak];
      func->SetParameter(pIndex, map.p0[ipeak]);
      func->SetParLimits(pIndex, fitLow, fitHigh);
      const PhotoPeakIndexedParConfig *entry =
        PhotoPeakIndexedEntry(config.multiP, ipeak);
      if (entry && entry->hasLimit) {
        func->SetParLimits(pIndex, entry->limitLow, entry->limitHigh);
      }
    }
  }

  if (map.relativeFwhm) {
    if (map.wScaleIndex >= 0) {
      func->SetParameter(map.wScaleIndex,
                         config.hasWScaleInit ? config.wScaleInit : 1.0);
      func->SetParLimits(map.wScaleIndex, 1.0e-6, 10.0);
      if (config.hasWScaleLimit) {
        func->SetParLimits(map.wScaleIndex,
                           config.wScaleLimitLow, config.wScaleLimitHigh);
      }
      if (config.hasWScaleFix) {
        func->FixParameter(map.wScaleIndex, config.wScaleFix);
      }
    }
  } else {
    for (int ipeak = 0; ipeak < map.nPeaks; ++ipeak) {
      if (map.hasWFixed[ipeak]) {
        continue;
      }
      const int wIndex = map.wIndex[ipeak];
      func->SetParameter(wIndex, map.w0[ipeak]);
      func->SetParLimits(wIndex, 1.0e-6, range);
      const PhotoPeakIndexedParConfig *entry =
        PhotoPeakIndexedEntry(config.multiW, ipeak);
      if (entry && entry->hasLimit) {
        func->SetParLimits(wIndex, entry->limitLow, entry->limitHigh);
      }
    }
  }

  for (int ipeak = 0; ipeak < map.nPeaks; ++ipeak) {
    const int hIndex = map.hIndex[ipeak];
    double hInit = h0[ipeak];
    const PhotoPeakIndexedParConfig *entry =
      PhotoPeakIndexedEntry(config.multiH, ipeak);
    if (entry && entry->hasInit) {
      hInit = entry->init;
    }
    func->SetParameter(hIndex, hInit);
    func->SetParLimits(hIndex, 0.0, hUpper);
    if (entry && entry->hasLimit) {
      func->SetParLimits(hIndex, entry->limitLow, entry->limitHigh);
    }
    if (entry && entry->hasFix) {
      func->FixParameter(hIndex, entry->fix);
    }
  }
  func->SetNpx(2000);
}

// ============== PhotoPeakGetPeakPosition ==============
// Purpose: Return a fitted or derived multi-peak centroid with uncertainty.
// Inputs: Peak index, TF1 pointer, and parameter map.
// Outputs: Derived value containing value, uncertainty, and source.
PhotoPeakDerivedValue PhotoPeakGetPeakPosition(int ipeak, TF1 *func,
                                               const PhotoPeakMultiParMap &map)
{
  PhotoPeakDerivedValue result;
  std::vector<double> par(map.parNames.size(), 0.0);
  for (int ipar = 0; ipar < static_cast<int>(par.size()); ++ipar) {
    par[ipar] = func->GetParameter(ipar);
  }
  result.value = PhotoPeakMultiPositionValue(ipeak, &par[0]);
  if (map.hasPFixed[ipeak] ||
      (map.relativePosition && map.hasFixedPositionAnchor)) {
    result.error = 0.0;
    result.source = "fixed";
  } else if (map.relativePosition) {
    result.error = func->GetParError(map.pIndex[0]);
    result.source = "derived from PMASTER";
  } else {
    result.error = func->GetParError(map.pIndex[ipeak]);
    result.source = "fit";
  }

  return result;
}

// ============== PhotoPeakGetPeakFwhm ==============
// Purpose: Return a fitted, fixed, or scaled multi-peak FWHM with uncertainty.
// Inputs: Peak index, TF1 pointer, and parameter map.
// Outputs: Derived value containing value, uncertainty, and source.
PhotoPeakDerivedValue PhotoPeakGetPeakFwhm(int ipeak, TF1 *func,
                                           const PhotoPeakMultiParMap &map)
{
  PhotoPeakDerivedValue result;
  std::vector<double> par(map.parNames.size(), 0.0);
  for (int ipar = 0; ipar < static_cast<int>(par.size()); ++ipar) {
    par[ipar] = func->GetParameter(ipar);
  }
  result.value = PhotoPeakMultiFwhmValue(ipeak, &par[0]);
  if (map.hasWFixed[ipeak]) {
    result.error = 0.0;
    result.source = "fixed";
  } else if (map.relativeFwhm) {
    if (map.wScaleIndex >= 0) {
      result.error = map.w0[ipeak] * func->GetParError(map.wScaleIndex);
      result.source = "derived from WSCALE";
    } else {
      result.error = 0.0;
      result.source = "fixed";
    }
  } else {
    result.error = func->GetParError(map.wIndex[ipeak]);
    result.source = "fit";
  }

  return result;
}

// ============== PhotoPeakMultiArea ==============
// Purpose: Calculate one multi-peak area from shared shape and peak values.
// Inputs: R, BETA, FWHM, height, and histogram bin width.
// Outputs: Photopeak area in histogram-bin counts.
double PhotoPeakMultiArea(double rPercent, double betaValue,
                          double wFwhm, double h, double binWidth)
{
  double par[kNPars] = {0.0};
  par[kR] = rPercent;
  par[kBeta] = betaValue;
  par[kW] = wFwhm;
  par[kH] = h;

  return PhotoPeakArea(par, binWidth);
}

// ============== PhotoPeakAreaDerivatives ==============
// Purpose: Calculate area derivatives with respect to R, BETA, W, and H.
// Inputs: R, BETA, FWHM, height, bin width, and derivative references.
// Outputs: Filled derivative values.
void PhotoPeakAreaDerivatives(double rPercent, double betaValue,
                              double wFwhmValue, double hValue,
                              double binWidth, double &dAdR,
                              double &dAdBeta, double &dAdW,
                              double &dAdH)
{
  const double invBinWidth = 1.0 / std::max(binWidth, 1.0e-12);
  const double r = rPercent / 100.0;
  const double beta = std::max(betaValue, 1.0e-12);
  const double wFwhm = std::max(wFwhmValue, 1.0e-12);
  const double h = hValue;
  const double y = wFwhm / (beta * 3.33021838);
  const double erfcY = TMath::Erfc(y);
  const double d = erfcY > 0.0 ? std::exp(-y * y) / erfcY : 0.0;
  dAdH = (r * 2.0 * beta * d + (1.0 - r) * wFwhm * 1.06446705) * invBinWidth;
  dAdR = 0.01 * h * (2.0 * beta * d - wFwhm * 1.06446705) * invBinWidth;
  dAdBeta = h * r * 2.0 * d *
    (1.0 + 2.0 * y * y - d * 1.12837917 * y) * invBinWidth;
  dAdW = h *
    (r * 2.0 * 0.600561216 * d * (d / 1.77245385 - y) +
     (1.0 - r) * 1.06446705) * invBinWidth;
}

// ============== PhotoPeakMultiAreaUncertainty ==============
// Purpose: Propagate full covariance to one multi-peak area uncertainty.
// Inputs: Peak index, TF1, covariance matrix, parameter map, and bin width.
// Outputs: One-sigma area uncertainty.
double PhotoPeakMultiAreaUncertainty(int ipeak, TF1 *func,
                                     const TMatrixDSym &cov,
                                     const PhotoPeakMultiParMap &map,
                                     double binWidth)
{
  std::vector<double> par(map.parNames.size(), 0.0);
  std::vector<double> grad(map.parNames.size(), 0.0);
  for (int ipar = 0; ipar < static_cast<int>(par.size()); ++ipar) {
    par[ipar] = func->GetParameter(ipar);
  }
  const double wFwhm = PhotoPeakMultiFwhmValue(ipeak, &par[0]);
  const double h = par[map.hIndex[ipeak]];
  double dAdR = 0.0;
  double dAdBeta = 0.0;
  double dAdW = 0.0;
  double dAdH = 0.0;
  PhotoPeakAreaDerivatives(par[kR], par[kBeta], wFwhm, h, binWidth,
                           dAdR, dAdBeta, dAdW, dAdH);
  if (map.useTail && !(par[kR] == 0.0)) {
    grad[kR] = dAdR;
    grad[kBeta] = dAdBeta;
  }
  if (map.hasWFixed[ipeak]) {
    dAdW = 0.0;
  } else if (map.relativeFwhm) {
    if (map.wScaleIndex >= 0) {
      grad[map.wScaleIndex] = dAdW * map.w0[ipeak];
    }
  } else {
    grad[map.wIndex[ipeak]] = dAdW;
  }
  grad[map.hIndex[ipeak]] = dAdH;

  double variance = 0.0;
  for (int i = 0; i < static_cast<int>(grad.size()); ++i) {
    for (int j = 0; j < static_cast<int>(grad.size()); ++j) {
      variance += grad[i] * cov(i, j) * grad[j];
    }
  }

  return variance > 0.0 ? std::sqrt(variance) : 0.0;
}

// ============== PhotoPeakPrintMultiResult ==============
// Purpose: Print grouped multi-peak values and raw fit parameters.
// Inputs: Histogram, fit metadata, TF1, covariance, map, and config.
// Outputs: Fit report printed to stdout.
void PhotoPeakPrintMultiResult(TH1 *hist, double fitLow, double fitHigh,
                               const std::vector<double> &peaks,
                               int mode,
                               const PhotoPeakTspectrumBgResult &bgResult,
                               const PhotoPeakFitTrial &trial,
                               TF1 *total,
                               const TMatrixDSym &cov,
                               const PhotoPeakMultiParMap &map,
                               const PhotoPeakFitConfig &config,
                               int nFitBins)
{
  std::printf("\nmultipeakfit result for %s\n", hist->GetName());
  std::printf("Fit range: [%g, %g], initial peaks:", fitLow, fitHigh);
  for (int ipeak = 0; ipeak < static_cast<int>(peaks.size()); ++ipeak) {
    std::printf(" %g", peaks[ipeak]);
  }
  std::printf("\n");
  std::printf("Requested mode: %s\n", PhotoPeakModeName(mode));
  std::printf("ROOT fit option: %s\n", config.rootFitOption.c_str());
  std::printf("TSpectrum background option: %s\n",
              PhotoPeakTspectrumBgModeName(bgResult.mode));
  if (bgResult.mode != kPhotoPeakTspectrumBgNone) {
    std::printf("TSpectrum iteration: %d\n", bgResult.iteration);
    std::printf("TSpectrum range: [%g, %g]\n",
                bgResult.rangeLow, bgResult.rangeHigh);
  }
  std::printf("Relative position fixed: %s\n",
              PhotoPeakBoolName(config.relativePosition));
  std::printf("Relative FWHM fixed: %s\n",
              PhotoPeakBoolName(config.relativeFwhm));
  std::printf("Fitting function: %s\n", trial.candidate.name);
  std::printf("Fit bins = %d, free parameters = %d\n", nFitBins, trial.freePars);
  std::printf("Fit status: %d\n", trial.status);
  if (trial.ndf > 0) {
    std::printf("chi2 = %.10g, ndf = %d, reduced chisq = %.10g\n",
                trial.chi2, trial.ndf, trial.reducedChi2);
  } else {
    std::printf("chi2 = %.10g, ndf = %d, reduced chisq = n/a\n",
                trial.chi2, trial.ndf);
  }

  std::printf("\nPeaks:\n");
  for (int ipeak = 0; ipeak < map.nPeaks; ++ipeak) {
    const PhotoPeakDerivedValue p = PhotoPeakGetPeakPosition(ipeak, total, map);
    const PhotoPeakDerivedValue w = PhotoPeakGetPeakFwhm(ipeak, total, map);
    const double h = total->GetParameter(map.hIndex[ipeak]);
    const double hErr = total->GetParError(map.hIndex[ipeak]);
    const double binWidth =
      hist->GetXaxis()->GetBinWidth(hist->GetXaxis()->FindFixBin(p.value));
    const double area =
      PhotoPeakMultiArea(total->GetParameter(kR), total->GetParameter(kBeta),
                         w.value, h, binWidth);
    const double areaErr =
      PhotoPeakMultiAreaUncertainty(ipeak, total, cov, map, binWidth);
    std::printf("  Peak %d:\n", ipeak);
    std::printf("    Position = % .10g +/- %.10g (%s)\n",
                p.value, p.error, p.source);
    std::printf("    Height   = % .10g +/- %.10g\n", h, hErr);
    std::printf("    FWHM     = % .10g +/- %.10g (%s)\n",
                w.value, w.error, w.source);
    std::printf("    Area     = % .10g +/- %.10g\n", area, areaErr);
  }

  std::printf("\nFit parameters:\n");
  for (int ipar = 0; ipar < total->GetNpar(); ++ipar) {
    std::printf("  %-8s = % .10g +/- %.10g\n",
                total->GetParName(ipar),
                total->GetParameter(ipar), total->GetParError(ipar));
  }
  std::printf("\n");
}

// ============== PhotoPeakCopyMultiParameters ==============
// Purpose: Copy fitted parameters from one TF1 to another.
// Inputs: Source TF1, destination TF1, and parameter count.
// Outputs: Destination TF1 parameter values.
void PhotoPeakCopyMultiParameters(TF1 *source, TF1 *dest, int npars)
{
  for (int ipar = 0; ipar < npars; ++ipar) {
    dest->SetParameter(ipar, source->GetParameter(ipar));
  }
}

// ============== PhotoPeakDrawMultiFit ==============
// Purpose: Draw multi-peak total, background, TSpectrum background, peaks, and legend.
// Inputs: Histogram, fit range, fitted function, background metadata, and map.
// Outputs: Updated canvas with fit curves.
void PhotoPeakDrawMultiFit(TH1 *hist, double fitLow, double fitHigh,
                           TF1 *total, const PhotoPeakTspectrumBgResult &bgResult,
                           const PhotoPeakMultiParMap &map)
{
  const int npars = static_cast<int>(map.parNames.size());
  if (!gPad) {
    new TCanvas("MultiPeak_canvas", "MultiPeak_canvas", 900, 650);
    hist->Draw();
  } else {
    gPad->cd();
    if (!hist->TestBit(TH1::kIsZoomed)) {
      hist->Draw();
    }
  }

  TF1 *totalDraw = total;
  if (bgResult.mode != kPhotoPeakTspectrumBgNone) {
    totalDraw = new TF1("MultiPeak_total_fit_display",
                        PhotoPeakMultiDisplayEval, fitLow, fitHigh, npars);
    PhotoPeakCopyMultiParameters(total, totalDraw, npars);
  }
  totalDraw->SetLineColor(kRed);
  totalDraw->SetLineStyle(1);
  totalDraw->SetLineWidth(4);
  totalDraw->SetNpx(2000);

  TF1 *bg = new TF1("MultiPeak_background_fit",
                    bgResult.mode == kPhotoPeakTspectrumBgNone ?
                    PhotoPeakMultiBgEval : PhotoPeakMultiBgDisplayEval,
                    fitLow, fitHigh, npars);
  PhotoPeakCopyMultiParameters(total, bg, npars);
  bg->SetLineColor(kBlack);
  bg->SetLineStyle(9);
  bg->SetLineWidth(4);
  bg->SetNpx(2000);

  TF1 *tspectrumBg = 0;
  if (bgResult.mode != kPhotoPeakTspectrumBgNone) {
    tspectrumBg = new TF1("MultiPeak_tspectrum_background_draw",
                          PhotoPeakTspectrumBgEval,
                          bgResult.rangeLow, bgResult.rangeHigh, 0);
    tspectrumBg->SetLineColor(kRed);
    tspectrumBg->SetLineStyle(9);
    tspectrumBg->SetLineWidth(3);
    tspectrumBg->SetNpx(2000);
  }

  const int colors[] = {
    kBlue + 1, kGreen + 2, kMagenta + 1, kCyan + 2,
    kOrange + 7, kViolet + 1, kAzure + 7, kSpring + 5,
    kPink + 7, kTeal + 3
  };
  TLegend *legend = new TLegend(0.12, 0.70, 0.42, 0.90);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.035);

  totalDraw->Draw("same");
  bg->Draw("same");
  if (tspectrumBg) {
    tspectrumBg->Draw("same");
  }

  for (int ipeak = 0; ipeak < map.nPeaks; ++ipeak) {
    TString name;
    name.Form("MultiPeak_peak_%d", ipeak);
    TF1 *peak = new TF1(name.Data(), PhotoPeakMultiPeakEval,
                        fitLow, fitHigh, npars + 1);
    PhotoPeakCopyMultiParameters(total, peak, npars);
    peak->SetParameter(npars, ipeak);
    peak->SetLineColor(colors[ipeak % 10]);
    peak->SetLineStyle(1);
    peak->SetLineWidth(2);
    peak->SetNpx(2000);
    peak->Draw("same");

    const PhotoPeakDerivedValue p =
      PhotoPeakGetPeakPosition(ipeak, total, map);
    TString label;
    label.Form("Peak %d: %.4g", ipeak, p.value);
    legend->AddEntry(peak, label.Data(), "l");
  }
  legend->Draw();
  if (gPad) {
    gPad->Modified();
    gPad->Update();
  }
}

// ============== PhotoPeakMultiFitRun ==============
// Purpose: Fit multiple photopeaks with shared GF3-style shape/background parameters.
// Inputs: Histogram pointer, fit range, initial peak x-values, and config.
// Outputs: Fit status; prints fit results and draws total/background/peak functions.
int PhotoPeakMultiFitRun(TH1 *hist, double fitLow, double fitHigh,
                         const std::vector<double> &inputPeaks,
                         const PhotoPeakFitConfig &config)
{
  if (!hist) {
    std::printf("multipeakfit ERROR: null histogram pointer.\n");
    return 1;
  }
  if (inputPeaks.empty()) {
    std::printf("multipeakfit ERROR: no input peaks.\n");
    return 2;
  }
  if (fitLow == fitHigh) {
    std::printf("multipeakfit ERROR: fitLow and fitHigh are equal.\n");
    return 3;
  }
  if (fitLow > fitHigh) {
    std::swap(fitLow, fitHigh);
  }

  std::vector<double> peaks = inputPeaks;
  for (int ipeak = 1; ipeak < static_cast<int>(peaks.size()); ++ipeak) {
    if (peaks[ipeak] < peaks[ipeak - 1]) {
      std::printf("multipeakfit WARNING: input peaks are not in increasing x order; preserving user order.\n");
      break;
    }
  }

  int mode = config.mode;
  if (mode != kPhotoPeakAuto &&
      mode != kPhotoPeakHighStat &&
      mode != kPhotoPeakLowStat) {
    std::printf("multipeakfit WARNING: unknown mode %d; using auto.\n", mode);
    mode = kPhotoPeakAuto;
  }

  PhotoPeakTspectrumBgResult bgResult =
    PhotoPeakPrepareTspectrumBackground(hist, fitLow, fitHigh, peaks[0], config);
  TH1 *fitHist = bgResult.fitHist;
  gPhotoPeakTspectrumBg = bgResult.mode == kPhotoPeakTspectrumBgNone ?
    0 : bgResult.bgHist;

  gPhotoPeakFitLow = fitLow;
  gPhotoPeakFitHigh = fitHigh;

  const double yLow = PhotoPeakBinContentAtX(fitHist, fitLow);
  const double yHigh = PhotoPeakBinContentAtX(fitHist, fitHigh);
  const double xMid = 0.5 * (fitLow + fitHigh);
  const double range = fitHigh - fitLow;
  const double a0 = 0.5 * (yLow + yHigh);
  const double b0 = (yHigh - yLow) / range;
  const double c0 = 0.0;
  const double r0 = 10.0;
  const double step0 = 0.25;
  const double maxInRange = PhotoPeakMaxInRange(fitHist, fitLow, fitHigh);
  const double hUpper = std::max(10.0 * maxInRange, 1.0);
  const int nFitBins = PhotoPeakBinCountInRange(fitHist, fitLow, fitHigh);

  std::vector<double> w0(peaks.size(), 0.0);
  std::vector<double> h0(peaks.size(), 0.0);
  for (int ipeak = 0; ipeak < static_cast<int>(peaks.size()); ++ipeak) {
    const PhotoPeakIndexedParConfig *wEntry =
      PhotoPeakIndexedEntry(config.multiW, ipeak);
    if (wEntry && wEntry->hasInit) {
      w0[ipeak] = wEntry->init;
    } else {
      w0[ipeak] = std::sqrt(std::max(9.0 + 0.004 * peaks[ipeak], 1.0e-12));
    }
    const double yPeak = PhotoPeakBinContentAtX(fitHist, peaks[ipeak]);
    const double linearBgAtPeak = a0 + b0 * (peaks[ipeak] - xMid);
    h0[ipeak] = std::max(yPeak - linearBgAtPeak, std::max(maxInRange, 1.0));
  }
  const double beta0 = 0.5 * w0[0];

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
  PhotoPeakMultiParMap maps[8];
  TMatrixDSym bestCov(1);
  TMatrixDSym fallbackCov(1);
  PhotoPeakMultiParMap bestMap;
  PhotoPeakMultiParMap fallbackMap;
  PhotoPeakFitTrial *bestTrial = 0;
  PhotoPeakFitTrial *fallbackTrial = 0;

  for (int i = 0; i < 8; ++i) {
    PhotoPeakFitTrial &trial = trials[i];
    maps[i] = PhotoPeakBuildMultiParMap(peaks, w0, trial.candidate, config);
    trial.freePars = PhotoPeakMultiFreeParameterCount(maps[i], trial.candidate, config);
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

    gPhotoPeakMultiMap = &maps[i];
    TString funcName;
    funcName.Form("MultiPeak_total_fit_%s", trial.candidate.name);
    trial.func = new TF1(funcName.Data(), PhotoPeakMultiEval,
                         fitLow, fitHigh,
                         static_cast<int>(maps[i].parNames.size()));
    PhotoPeakConfigureMultiFunction(trial.func, maps[i], trial.candidate,
                                    a0, b0, c0, r0, beta0, step0, h0,
                                    fitLow, fitHigh, range, hUpper, config);
    TFitResultPtr result = fitHist->Fit(trial.func, config.rootFitOption.c_str());
    trial.status = int(result);
    trial.chi2 = trial.func->GetChisquare();
    trial.ndf = trial.func->GetNDF();
    trial.reducedChi2 = trial.ndf > 0 ? trial.chi2 / trial.ndf :
      std::numeric_limits<double>::quiet_NaN();
    TMatrixDSym currentCov(static_cast<int>(maps[i].parNames.size()));
    if (result.Get()) {
      currentCov = result->GetCovarianceMatrix();
    }
    if (isBasic) {
      fallbackTrial = &trial;
      fallbackCov.ResizeTo(currentCov);
      fallbackCov = currentCov;
      fallbackMap = maps[i];
    }
    if (mode == kPhotoPeakAuto) {
      if (PhotoPeakTrialIsBetter(&trial, bestTrial)) {
        bestTrial = &trial;
        bestCov.ResizeTo(currentCov);
        bestCov = currentCov;
        bestMap = maps[i];
      }
    } else {
      bestTrial = &trial;
      bestCov.ResizeTo(currentCov);
      bestCov = currentCov;
      bestMap = maps[i];
    }
  }

  if (!bestTrial) {
    bestTrial = fallbackTrial;
    bestCov.ResizeTo(fallbackCov);
    bestCov = fallbackCov;
    bestMap = fallbackMap;
  }
  if (!bestTrial || !bestTrial->func) {
    std::printf("multipeakfit ERROR: no fitting function was available.\n");
    return 4;
  }

  for (int i = 0; i < 8; ++i) {
    if (&trials[i] != bestTrial && trials[i].func) {
      delete trials[i].func;
      trials[i].func = 0;
    }
  }

  gPhotoPeakMultiMapStorage = bestMap;
  gPhotoPeakMultiMap = &gPhotoPeakMultiMapStorage;
  PhotoPeakPrintMultiResult(hist, fitLow, fitHigh, peaks, mode, bgResult,
                            *bestTrial, bestTrial->func, bestCov,
                            gPhotoPeakMultiMapStorage, config, nFitBins);
  PhotoPeakDrawMultiFit(hist, fitLow, fitHigh, bestTrial->func, bgResult,
                        gPhotoPeakMultiMapStorage);

  return bestTrial->status;
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
  if (config.tspectrumBgMode != kPhotoPeakTspectrumBgNone &&
      config.tspectrumBgMode != kPhotoPeakTspectrumBgGlobal &&
      config.tspectrumBgMode != kPhotoPeakTspectrumBgLocal) {
    std::printf("photopeakfit WARNING: unknown TSpectrum background option %d; using none.\n",
                config.tspectrumBgMode);
  }

  PhotoPeakTspectrumBgResult bgResult =
    PhotoPeakPrepareTspectrumBackground(hist, fitLow, fitHigh, peak0, config);
  TH1 *fitHist = bgResult.fitHist;
  gPhotoPeakTspectrumBg = bgResult.mode == kPhotoPeakTspectrumBgNone ?
    0 : bgResult.bgHist;

  gPhotoPeakFitLow = fitLow;
  gPhotoPeakFitHigh = fitHigh;

  const double yLow = PhotoPeakBinContentAtX(fitHist, fitLow);
  const double yHigh = PhotoPeakBinContentAtX(fitHist, fitHigh);
  const double yPeak = PhotoPeakBinContentAtX(fitHist, peak0);
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
  const double maxInRange = PhotoPeakMaxInRange(fitHist, fitLow, fitHigh);
  const double h0 = std::max(yPeak - linearBgAtPeak, std::max(maxInRange, 1.0));
  const double hUpper = std::max(10.0 * maxInRange, h0 * 10.0);
  const int nFitBins = PhotoPeakBinCountInRange(fitHist, fitLow, fitHigh);
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
    TFitResultPtr result = fitHist->Fit(trial.func, config.rootFitOption.c_str());
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
  std::printf("ROOT fit option: %s\n", config.rootFitOption.c_str());
  std::printf("TSpectrum background option: %s\n",
              PhotoPeakTspectrumBgModeName(bgResult.mode));
  if (bgResult.mode != kPhotoPeakTspectrumBgNone) {
    std::printf("TSpectrum iteration: %d\n", bgResult.iteration);
    std::printf("TSpectrum range: [%g, %g]\n",
                bgResult.rangeLow, bgResult.rangeHigh);
  }
  std::printf("Relative position fixed: false\n");
  std::printf("Relative FWHM fixed: false\n");
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

  TF1 *totalDraw = total;
  if (bgResult.mode != kPhotoPeakTspectrumBgNone) {
    totalDraw = new TF1("PhotoPeak_total_fit_display",
                        PhotoPeakDisplayEval, fitLow, fitHigh, kNPars);
    totalDraw->SetParameters(par);
  }
  totalDraw->SetLineColor(kRed);
  totalDraw->SetLineStyle(1);
  totalDraw->SetLineWidth(3);
  totalDraw->SetNpx(2000);

  TF1 *bg = new TF1("PhotoPeak_background_fit",
                    bgResult.mode == kPhotoPeakTspectrumBgNone ?
                    PhotoPeakBgEval : PhotoPeakBgDisplayEval,
                    fitLow, fitHigh, kNPars);
  bg->SetParameters(par);
  bg->SetLineColor(kBlack);
  bg->SetLineStyle(2);
  bg->SetLineWidth(3);
  bg->SetNpx(2000);

  TF1 *tspectrumBg = 0;
  if (bgResult.mode != kPhotoPeakTspectrumBgNone) {
    tspectrumBg = new TF1("PhotoPeak_tspectrum_background_draw",
                          PhotoPeakTspectrumBgEval,
                          bgResult.rangeLow, bgResult.rangeHigh, 0);
    tspectrumBg->SetLineColor(kRed);
    tspectrumBg->SetLineStyle(7);
    tspectrumBg->SetLineWidth(3);
    tspectrumBg->SetNpx(2000);
  }

  if (!gPad) {
    new TCanvas("PhotoPeak_canvas", "PhotoPeak_canvas", 900, 650);
    hist->Draw();
  } else {
    gPad->cd();
    if (!hist->TestBit(TH1::kIsZoomed)) {
      hist->Draw();
    }
  }
  totalDraw->Draw("same");
  bg->Draw("same");
  if (tspectrumBg) {
    tspectrumBg->Draw("same");
  }
  if (gPad) {
    gPad->Modified();
    gPad->Update();
  }

  return status;
}

// ============== PhotoPeakConfigPeakPosition ==============
// Purpose: Read the single-peak initial position from a config object.
// Inputs: Config object and config file label.
// Outputs: Peak position through reference; zero on success.
int PhotoPeakConfigPeakPosition(const PhotoPeakFitConfig &config,
                                const char *configFile,
                                double &peak0)
{
  if (!config.hasInit[kP]) {
    std::printf("photopeakfit ERROR: config file '%s' must define init P when peak0 is omitted.\n",
                configFile ? configFile : "");
    return 1;
  }
  peak0 = config.init[kP];

  return 0;
}

// ============== PhotoPeakConfigPeakList ==============
// Purpose: Read a contiguous multi-peak position list from a config object.
// Inputs: Config object, config file label, and output peak vector.
// Outputs: Zero on success.
int PhotoPeakConfigPeakList(const PhotoPeakFitConfig &config,
                            const char *configFile,
                            std::vector<double> &peaks)
{
  peaks.clear();
  if (config.multiP.empty() || !config.multiP[0].hasInit) {
    std::printf("multipeakfit ERROR: config file '%s' must define init P[0] when peaks are omitted.\n",
                configFile ? configFile : "");
    return 1;
  }

  bool foundGap = false;
  int gapIndex = -1;
  for (int ipeak = 0; ipeak < static_cast<int>(config.multiP.size()); ++ipeak) {
    if (config.multiP[ipeak].hasInit) {
      if (foundGap) {
        std::printf("multipeakfit ERROR: config file '%s' defines init P[%d] but is missing init P[%d].\n",
                    configFile ? configFile : "", ipeak, gapIndex);
        return 2;
      }
      peaks.push_back(config.multiP[ipeak].init);
    } else if (!foundGap) {
      foundGap = true;
      gapIndex = ipeak;
    }
  }
  if (peaks.empty()) {
    std::printf("multipeakfit ERROR: config file '%s' does not define any init P[i] values.\n",
                configFile ? configFile : "");
    return 3;
  }

  return 0;
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

// ============== photopeakfit ==============
// Purpose: Fit one photopeak using a config-defined initial peak position.
// Inputs: Histogram pointer, lower x range, upper x range, and config file.
// Outputs: Fit status; prints fit results and draws total/background functions.
int photopeakfit(TH1 *hist, double fitLow, double fitHigh,
                 const char *configFile)
{
  PhotoPeakFitConfig config;
  const int readStatus = PhotoPeakReadConfig(configFile, config);
  if (readStatus != 0) {
    return 10 + readStatus;
  }
  double peak0 = 0.0;
  const int peakStatus =
    PhotoPeakConfigPeakPosition(config, configFile, peak0);
  if (peakStatus != 0) {
    return 20 + peakStatus;
  }

  return PhotoPeakFitRun(hist, fitLow, fitHigh, peak0, config);
}

// ============== multipeakfit ==============
// Purpose: Fit multiple photopeaks with default multi-peak configuration.
// Inputs: Histogram pointer, fit range, and initial peak x-values.
// Outputs: Fit status; prints fit results and draws total/background/peak functions.
int multipeakfit(TH1 *hist, double fitLow, double fitHigh,
                 std::initializer_list<double> peaks)
{
  PhotoPeakFitConfig config;
  std::vector<double> peakVector(peaks.begin(), peaks.end());

  return PhotoPeakMultiFitRun(hist, fitLow, fitHigh, peakVector, config);
}

// ============== multipeakfit ==============
// Purpose: Fit multiple photopeaks using a text configuration file.
// Inputs: Histogram pointer, fit range, initial peak x-values, and config file.
// Outputs: Fit status; prints fit results and draws total/background/peak functions.
int multipeakfit(TH1 *hist, double fitLow, double fitHigh,
                 std::initializer_list<double> peaks,
                 const char *configFile)
{
  PhotoPeakFitConfig config;
  const int readStatus = PhotoPeakReadConfig(configFile, config);
  if (readStatus != 0) {
    return 10 + readStatus;
  }
  std::vector<double> peakVector(peaks.begin(), peaks.end());

  return PhotoPeakMultiFitRun(hist, fitLow, fitHigh, peakVector, config);
}

// ============== multipeakfit ==============
// Purpose: Fit multiple photopeaks using config-defined initial peak positions.
// Inputs: Histogram pointer, fit range, and config file.
// Outputs: Fit status; prints fit results and draws total/background/peak functions.
int multipeakfit(TH1 *hist, double fitLow, double fitHigh,
                 const char *configFile)
{
  PhotoPeakFitConfig config;
  const int readStatus = PhotoPeakReadConfig(configFile, config);
  if (readStatus != 0) {
    return 10 + readStatus;
  }
  std::vector<double> peakVector;
  const int peakStatus =
    PhotoPeakConfigPeakList(config, configFile, peakVector);
  if (peakStatus != 0) {
    return 20 + peakStatus;
  }

  return PhotoPeakMultiFitRun(hist, fitLow, fitHigh, peakVector, config);
}

// ============== multipeakfit ==============
// Purpose: Fit multiple photopeaks with a C++ configuration object.
// Inputs: Histogram pointer, fit range, initial peak x-values, and config.
// Outputs: Fit status; prints fit results and draws total/background/peak functions.
int multipeakfit(TH1 *hist, double fitLow, double fitHigh,
                 const std::vector<double> &peaks,
                 const PhotoPeakFitConfig &config)
{
  return PhotoPeakMultiFitRun(hist, fitLow, fitHigh, peaks, config);
}

// ============== multipeakfit ==============
// Purpose: Fit multiple photopeaks with a C-style peak array and configuration.
// Inputs: Histogram pointer, fit range, peak count, peak x-array, and config.
// Outputs: Fit status; prints fit results and draws total/background/peak functions.
int multipeakfit(TH1 *hist, double fitLow, double fitHigh,
                 int nPeaks, const double *peaks,
                 const PhotoPeakFitConfig &config)
{
  if (nPeaks <= 0 || !peaks) {
    std::printf("multipeakfit ERROR: invalid peak array.\n");
    return 2;
  }
  std::vector<double> peakVector(peaks, peaks + nPeaks);

  return PhotoPeakMultiFitRun(hist, fitLow, fitHigh, peakVector, config);
}
