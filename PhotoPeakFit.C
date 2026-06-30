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

enum PhotoPeakPar {
  kA = 0,
  kB = 1,
  kC = 2,
  kR = 3,
  kBeta = 4,
  kStep = 5,
  kP = 6,
  kW = 7,
  kH = 8,
  kNPars = 9
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
// Purpose: Fit one histogram photopeak with a RadWare/GF3-style function.
// Inputs: Histogram pointer, lower x range, upper x range, initial peak x.
// Outputs: Fit status; prints fit results and draws total/background functions.
int photopeakfit(TH1 *hist, double fitLow, double fitHigh, double peak0)
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

  TF1 *total = new TF1("PhotoPeak_total_fit", PhotoPeakEval, fitLow, fitHigh, kNPars);
  total->SetParNames("bg0", "bg1", "bg2", "R", "BETA", "STEP",
                     "Centroid", "FWHM", "Height");
  total->SetParameters(a0, b0, c0, r0, beta0, step0, peak0, w0, h0);
  total->SetParLimits(kR, 0.0, 100.0);
  total->SetParLimits(kBeta, 1.0e-6, 10.0 * range);
  total->SetParLimits(kStep, 0.0, 100.0);
  total->SetParLimits(kP, fitLow, fitHigh);
  total->SetParLimits(kW, 1.0e-6, range);
  total->SetParLimits(kH, 0.0, hUpper);
  total->SetLineColor(kRed);
  total->SetLineStyle(1);
  total->SetLineWidth(2);
  total->SetNpx(2000);

  TFitResultPtr result = hist->Fit(total, "RQS");
  const int status = int(result);
  const double chi2 = total->GetChisquare();
  const int ndf = total->GetNDF();
  const double reducedChi2 = ndf > 0 ? chi2 / ndf : 0.0;

  double par[kNPars];
  for (int i = 0; i < kNPars; ++i) {
    par[i] = total->GetParameter(i);
  }

  TMatrixDSym cov(kNPars);
  if (result.Get()) {
    cov = result->GetCovarianceMatrix();
  }
  const double binWidth = hist->GetXaxis()->GetBinWidth(hist->GetXaxis()->FindFixBin(par[kP]));
  const double area = PhotoPeakArea(par, binWidth);
  const double areaErr = PhotoPeakAreaUncertainty(par, cov, binWidth);
  const int printOrder[kNPars] = {kP, kH, kW, kR, kBeta, kStep, kA, kB, kC};

  std::printf("\nphotopeakfit result for %s\n", hist->GetName());
  std::printf("Fit range: [%g, %g], initial peak position: %g\n", fitLow, fitHigh, peak0);
  std::printf("Fit status: %d\n", status);
  std::printf("chi2 = %.10g, ndf = %d, reduced chisq = %.10g\n", chi2, ndf, reducedChi2);
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
