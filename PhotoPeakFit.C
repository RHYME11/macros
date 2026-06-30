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
#include <cmath>
#include <cstdio>
#include <limits>

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

double gPhotoPeakFitLow = 0.0;
double gPhotoPeakFitHigh = 0.0;

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
    return "high_stat";
  }
  if (mode == kPhotoPeakLowStat) {
    return "low_stat";
  }

  return "auto";
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
                                double range, double hUpper)
{
  func->SetParNames("bg0", "bg1", "bg2", "R", "BETA", "STEP",
                    "Centroid", "FWHM", "Height");
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
  func->SetNpx(2000);
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

// ============== photopeakfit ==============
// Purpose: Fit one histogram photopeak with a selected RadWare/GF3-style function.
// Inputs: Histogram pointer, lower x range, upper x range, initial peak x, and mode.
// Outputs: Fit status; prints fit results and draws total/background functions.
int photopeakfit(TH1 *hist, double fitLow, double fitHigh, double peak0,
                 int mode = kPhotoPeakAuto)
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
    trial.freePars = PhotoPeakFreeParameterCount(trial.candidate);
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
                               range, hUpper);
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
