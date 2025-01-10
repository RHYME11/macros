
// ==== Polynomial ===== //
double Poly(double *dim, double *par, int order){
  double result = 0.0;
  for(int i=0;i<=order;i++){
    result += *(par+i) *TMath::Power(dim[0],i); 
  }
  return result;
}

// ==== StepFunction ==== //
double StepFunction(double *dim, double *par){
  double x = dim[0];

  double height = par[0];
  double cent   = par[1];
  double sigma  = par[2];
  double step   = par[3];

  return height*(step/100.0) *TMath::Erfc((x-cent)/(TMath::Sqrt(2.0)*sigma));
}

// ==== Step BG ==== //
double StepBG(double *dim, double *par){
  return StepFunction(dim,par) + Poly(dim, (par+4), 0);
}


// ==== multiple Gaus ==== //
// gaus_order = 1 ==> single Gaus;
// gaus_order = 2 ==> double Gaus;
double MultGaus(double *dim, double *par, int gaus_order){
   
  double result = 0.0;
  double x      = dim[0];
  double sigma  = par[0];
  for(int i=0;i<gaus_order;i++){
    double height = par[1+2*i];
    double cent   = par[2+2*i];
    result += height *TMath::Gaus(x,cent,sigma);
  }  
  
  return result;
}


// ==== single Gaus ==== //
double Gaus(double *dim, double *par) {
  // - dim[0]: channels to fit
  // - par[0]: height of peak
  // - par[1]: cent of peak
  // - par[2]: sigma
  // - par[3]: relative height of skewed gaus to gaus

  double x      = dim[0];
  double height = par[0];
  double cent   = par[1];
  double sigma  = par[2];
  double R      = par[3];

  return height*(1.0-R/100.0)*TMath::Gaus(x,cent,sigma);
}



// ==== SkewedGaus ==== //
double SkewedGaus(double *dim, double*par){
  // StepFunction(dim,par) + PolyBg
  // - par[0]: height of peak
  // - par[1]: cent of peak
  // - par[2]: sigma
  // - par[3]: relative height of skewed gaus to gaus
  // - par[4]: "skewedness" of the skewed gaussin

  Double_t x      = dim[0]; //channel number used for fitting
  Double_t height = par[0]; //height of photopeak
  Double_t cent   = par[1]; //Peak Centroid of non skew gaus
  Double_t sigma  = par[2]; //standard deviation of gaussian
  Double_t R      = par[3]; //relative height of skewed gaussian
  Double_t beta   = par[4]; //"skewedness" of the skewed gaussian

  double scaling = R*height/100.0;

  double fterm = (x-cent)/(sigma*TMath::Sqrt(2.));
  double sterm = sigma /  (beta *TMath::Sqrt(2.));

  return scaling * TMath::Exp((x-cent)/beta) * TMath::Erfc(fterm + sterm); 
}


// ==== PhotoPeak ==== //
double PhotoPeak(double *dim, double *par){
  return Gaus(dim,par) + SkewedGaus(dim,par);
}


// ==== PhotoPeak + StepBG ==== //
double PhotoPeakBG(double *dim, double *par){
  // - dim[0]: channels to fit
  // - par[0]: height of peak
  // - par[1]: cent of peak
  // - par[2]: sigma
  // - par[3]: relative height of skewed gaus to gaus
  // - par[4]: "skewedness" of the skewed gaussin
  // - par[5]: size of stepfunction step.

  // - par[6]: base bg height.
  // - par[7]: slope of bg.
  
  double spar[4];
  spar[0] = par[0];
  spar[1] = par[1];
  spar[2] = par[2];
  spar[3] = par[5];  //stepsize;
  return Gaus(dim,par) + SkewedGaus(dim,par) + StepFunction(dim,spar) + Poly(dim,par+6,0);
}


// ==== PhotoPeak + StepBG if out of range ==== //
double PhotoPeakBGExcludeRegion(double *dim, double *par){
  // - dim[0]: channels to fit
  // - par[0]: height of peak
  // - par[1]: cent of peak
  // - par[2]: sigma
  // - par[3]: relative height of skewed gaus to gaus
  // - par[4]: "skewedness" of the skewed gaussin
  // - par[5]: size of stepfunction step.
  
  // - par[6]: base bg height.
  
  // - par[7]: exclude low;
  // - par[8]: exclude high;

  if(dim[0]>par[7] && dim[0]<par[8]) {
    TF1::RejectPoint();
    return 0;
  }
  double spar[4];
  spar[0] = par[0];
  spar[1] = par[1];
  spar[2] = par[2];
  spar[3] = par[5];  //stepsize;
  return Gaus(dim,par) + SkewedGaus(dim,par) + StepFunction(dim,spar) + Poly(dim,par+6,0);
}


//==============================================================================================//
//All Bg is initialized as Poly2;
//Fix the last parameter = 0 ==> obtain linear bg 
double SingleGausBG(double *dim, double *par){
  double result = MultGaus(dim, par, 1);
  double bg     = Poly(dim, par+3, 2);
  
  return result + bg;  
}

double DoubleGausBG(double *dim, double *par){
  double result = MultGaus(dim, par, 2);
  double bg     = Poly(dim, par+5, 2);
  
  return result + bg;  
}

double TripleGausBG(double *dim, double *par){
  double result = MultGaus(dim, par, 3);
  double bg     = Poly(dim, par+7, 2);
  
  return result + bg;  
}

double QuadGausBG(double *dim, double *par){
  double result = MultGaus(dim, par, 4);
  double bg     = Poly(dim, par+9, 2);
  
  return result + bg;  
}

//==============================================================================================//




// ==== Fitting Functions ==== //
void SingleGausFit(TH1 *hist, double lower, double upper, double bgquad=0.0){
  hist->GetListOfFunctions()->Clear();

  TF1 *fx = new TF1("fx", SingleGausBG, lower, upper, 6);
  TF1 *fbg= new TF1("fbg","pol2",lower, upper);

  fx->SetParName(0, "Sigma");
  fx->SetParName(1, "Height");
  fx->SetParName(2, "Centroid");
  fx->SetParName(3, "BG offset");
  fx->SetParName(4, "BG slope");
  fx->SetParName(5, "BG Quad");

  fx->SetParameter(0, 2);
  fx->SetParameter(1, hist->GetBinContent(hist->FindBin((lower+upper)/2)));
  fx->SetParameter(2, (lower+upper)/2.);
  fx->SetParameter(3, hist->GetBinContent(hist->FindBin(lower)));
  fx->SetParameter(4, -0.1);
  fx->SetParameter(5,bgquad);

  fx->SetParLimits(0, 0, 5);
  fx->SetParLimits(1, 0, (hist->GetMaximum()*10));
  fx->SetParLimits(3, -fabs(hist->GetMaximum()*10), fabs(hist->GetMaximum()*10));
  fx->SetParLimits(4,-100,100);
  fx->SetParLimits(5,-100,0);

  int npar = fx->GetNpar();
  if(bgquad == 0.0){
    fx->FixParameter((npar-1),0);
  }
  hist->Fit(fx,"","",lower,upper);
  
  double bgpar[3] = {fx->GetParameter(npar-3), fx->GetParameter(npar-2), fx->GetParameter(npar-1)};
  fbg->SetParameters(bgpar);
  fbg->SetLineColor(kBlack);
  fbg->SetLineStyle(9);

  double rebin = hist->GetBinWidth(hist->FindBin(lower));
  double bg = fbg->Integral(lower,upper)/rebin;
  double sum = hist->Integral(hist->FindBin(lower), hist->FindBin(upper)-1)-bg;
  double area = fx->Integral(lower,upper)/rebin - bg;
  double chi2 = fx->GetChisquare()/(fx->GetNDF()-1);
	hist->GetListOfFunctions()->Add(fbg);
  system("Color 02");
  std::cout<< "Name: " << "gausbg" << std::endl;
  std::cout<< "Area: " << area  << std::endl;
  std::cout<< "Sum:  " << sum  << std::endl;
  std::cout<< "Chi2: " << chi2 << std::endl;
  gPad->Modified();
  gPad->Update();
}

void DoubleGausFit(TH1 *hist, double cent1, double cent2, double lower, double upper, double bgquad=0.0){
  hist->GetListOfFunctions()->Clear();

  TF1 *fx = new TF1("fx", DoubleGausBG, lower, upper, 8);
  TF1 *f1 = new TF1("f1", SingleGausBG, lower, upper, 6);
  TF1 *f2 = new TF1("f2", SingleGausBG, lower, upper, 6);
  TF1 *fbg= new TF1("fbg","pol2",lower, upper);

  fx->SetParName(0, "Sigma");
  fx->SetParName(1, "Height1");
  fx->SetParName(2, "Centroid1");
  fx->SetParName(3, "Height2");
  fx->SetParName(4, "Centroid2");
  fx->SetParName(5, "BG offset");
  fx->SetParName(6, "BG slope");
  fx->SetParName(7, "BG Quad");

  fx->SetParameter(0, 2);
  fx->SetParameter(1, hist->GetBinContent(hist->FindBin(cent1)));
  fx->SetParameter(2, cent1);
  fx->SetParameter(3, hist->GetBinContent(hist->FindBin(cent2)));
  fx->SetParameter(4, cent2);
  fx->SetParameter(5, hist->GetBinContent(hist->FindBin(lower)));
  fx->SetParameter(6, -0.1);
  fx->SetParameter(7,bgquad);

  fx->SetParLimits(0, 0, 5);
  fx->SetParLimits(1, 0, (hist->GetMaximum()*10));
  fx->SetParLimits(2, lower, upper);
  fx->SetParLimits(3, 0, (hist->GetMaximum()*10));
  fx->SetParLimits(4, lower, upper);
  fx->SetParLimits(5, -fabs(hist->GetMaximum()*10), fabs(hist->GetMaximum()*10));
  fx->SetParLimits(6,-100,100);
  fx->SetParLimits(7,-100,0);

  int npar = fx->GetNpar();
  if(bgquad == 0.0){
    fx->FixParameter((npar-1),0);
  }
  hist->Fit(fx,"","",lower,upper);
  
  double bgpar[3] = {fx->GetParameter(npar-3), fx->GetParameter(npar-2), fx->GetParameter(npar-1)};
  fbg->SetParameters(bgpar);
  fbg->SetLineColor(kBlack);
  fbg->SetLineStyle(9);
  f1->SetParameter(0, fx->GetParameter(0));
  f1->SetParameter(1, fx->GetParameter(1));
  f1->SetParameter(2, fx->GetParameter(2));
  f1->SetParameter(3, bgpar[0]);
  f1->SetParameter(4, bgpar[1]);
  f1->SetParameter(5, bgpar[2]);
  f1->SetLineColor(kBlue);
  f2->SetParameter(0, fx->GetParameter(0));
  f2->SetParameter(1, fx->GetParameter(3));
  f2->SetParameter(2, fx->GetParameter(4));
  f2->SetParameter(3, bgpar[0]);
  f2->SetParameter(4, bgpar[1]);
  f2->SetParameter(5, bgpar[2]);
  f2->SetLineColor(kGreen);

  hist->GetListOfFunctions()->Add(f1);
  hist->GetListOfFunctions()->Add(f2);
  hist->GetListOfFunctions()->Add(fbg);

  double rebin = hist->GetBinWidth(hist->FindBin(lower));
  double bg    = fbg->Integral(lower,upper)/rebin;
  double sum   = hist->Integral(hist->FindBin(lower), hist->FindBin(upper)-1)-bg;
  double area  = fx->Integral(lower,upper)/rebin - bg;
  double area1 = f1->Integral(lower,upper)/rebin - bg;
  double area2 = f2->Integral(lower,upper)/rebin - bg;
  double chi2  = fx->GetChisquare()/(fx->GetNDF()-1);
	hist->GetListOfFunctions()->Add(fbg);
  system("Color 02");
  std::cout<< "Name: "  << "gausbg" << std::endl;
  std::cout<< "Area: "  << area  << std::endl;
  std::cout<< "Area1: " << area1 << std::endl;
  std::cout<< "Area2: " << area2 << std::endl;
  std::cout<< "Sum:  "  << sum  << std::endl;
  std::cout<< "Chi2: "  << chi2 << std::endl;
  gPad->Modified();
  gPad->Update();
}


void TripleGausFit(TH1 *hist, double cent1, double cent2, double cent3,
                              double lower, double upper, double bgquad=0.0){
  hist->GetListOfFunctions()->Clear();

  TF1 *fx = new TF1("fx", TripleGausBG, lower, upper, 10);
  TF1 *f1 = new TF1("f1", SingleGausBG, lower, upper, 6);
  TF1 *f2 = new TF1("f2", SingleGausBG, lower, upper, 6);
  TF1 *f3 = new TF1("f3", SingleGausBG, lower, upper, 6);
  TF1 *fbg= new TF1("fbg","pol2",lower, upper);

  fx->SetParName(0, "Sigma");
  fx->SetParName(1, "Height1");
  fx->SetParName(2, "Centroid1");
  fx->SetParName(3, "Height2");
  fx->SetParName(4, "Centroid2");
  fx->SetParName(5, "Height3");
  fx->SetParName(6, "Centroid3");
  fx->SetParName(7, "BG offset");
  fx->SetParName(8, "BG slope");
  fx->SetParName(9, "BG Quad");

  fx->SetParameter(0, 2);
  fx->SetParameter(1, hist->GetBinContent(hist->FindBin(cent1)));
  fx->SetParameter(2, cent1);
  fx->SetParameter(3, hist->GetBinContent(hist->FindBin(cent2)));
  fx->SetParameter(4, cent2);
  fx->SetParameter(5, hist->GetBinContent(hist->FindBin(cent3)));
  fx->SetParameter(6, cent3);
  fx->SetParameter(7, hist->GetBinContent(hist->FindBin(lower)));
  fx->SetParameter(8, -0.1);
  fx->SetParameter(9,bgquad);

  fx->SetParLimits(0, 0, 5);
  fx->SetParLimits(1, 0, (hist->GetMaximum()*10));
  fx->SetParLimits(3, 0, (hist->GetMaximum()*10));
  fx->SetParLimits(5, 0, (hist->GetMaximum()*10));
  fx->SetParLimits(7, -fabs(hist->GetMaximum()*10), fabs(hist->GetMaximum()*10));
  fx->SetParLimits(8,-100,100);
  fx->SetParLimits(9,-100,0);

  int npar = fx->GetNpar();
  if(bgquad == 0.0){
    fx->FixParameter((npar-1),0);
  }
  hist->Fit(fx,"","",lower,upper);
  
  double bgpar[3] = {fx->GetParameter(npar-3), fx->GetParameter(npar-2), fx->GetParameter(npar-1)};
  fbg->SetParameters(bgpar);
  fbg->SetLineColor(kBlack);
  fbg->SetLineStyle(9);
  f1->SetParameter(0, fx->GetParameter(0));
  f1->SetParameter(1, fx->GetParameter(1));
  f1->SetParameter(2, fx->GetParameter(2));
  f1->SetParameter(3, bgpar[0]);
  f1->SetParameter(4, bgpar[1]);
  f1->SetParameter(5, bgpar[2]);
  f1->SetLineColor(kBlue);
  f2->SetParameter(0, fx->GetParameter(0));
  f2->SetParameter(1, fx->GetParameter(3));
  f2->SetParameter(2, fx->GetParameter(4));
  f2->SetParameter(3, bgpar[0]);
  f2->SetParameter(4, bgpar[1]);
  f2->SetParameter(5, bgpar[2]);
  f2->SetLineColor(kGreen);
  f3->SetParameter(0, fx->GetParameter(0));
  f3->SetParameter(1, fx->GetParameter(5));
  f3->SetParameter(2, fx->GetParameter(6));
  f3->SetParameter(3, bgpar[0]);
  f3->SetParameter(4, bgpar[1]);
  f3->SetParameter(5, bgpar[2]);
  f3->SetLineColor(kMagenta);

  hist->GetListOfFunctions()->Add(f1);
  hist->GetListOfFunctions()->Add(f2);
  hist->GetListOfFunctions()->Add(f3);
  hist->GetListOfFunctions()->Add(fbg);

  double rebin = hist->GetBinWidth(hist->FindBin(lower));
  double bg    = fbg->Integral(lower,upper)/rebin;
  double sum   = hist->Integral(hist->FindBin(lower), hist->FindBin(upper)-1)-bg;
  double area  = fx->Integral(lower,upper)/rebin - bg;
  double area1 = f1->Integral(lower,upper)/rebin - bg;
  double area2 = f2->Integral(lower,upper)/rebin - bg;
  double area3 = f3->Integral(lower,upper)/rebin - bg;
  double chi2  = fx->GetChisquare()/(fx->GetNDF()-1);
	hist->GetListOfFunctions()->Add(fbg);
  system("Color 02");
  std::cout<< "Name: "  << "gausbg" << std::endl;
  std::cout<< "Area: "  << area  << std::endl;
  std::cout<< "Area1: " << area1 << std::endl;
  std::cout<< "Area2: " << area2 << std::endl;
  std::cout<< "Area3: " << area3 << std::endl;
  std::cout<< "Sum:  "  << sum  << std::endl;
  std::cout<< "Chi2: "  << chi2 << std::endl;
  gPad->Modified();
  gPad->Update();
}


void QuadGausFit(TH1 *hist, double cent1, double cent2, double cent3, double cent4,
                            double lower, double upper, double bgquad=0.0){
  hist->GetListOfFunctions()->Clear();

  TF1 *fx = new TF1("fx", QuadGausBG, lower, upper, 12);
  TF1 *f1 = new TF1("f1", SingleGausBG, lower, upper, 6);
  TF1 *f2 = new TF1("f2", SingleGausBG, lower, upper, 6);
  TF1 *f3 = new TF1("f3", SingleGausBG, lower, upper, 6);
  TF1 *f4 = new TF1("f4", SingleGausBG, lower, upper, 6);
  TF1 *fbg= new TF1("fbg","pol2",lower, upper);

  fx->SetParName(0, "Sigma");
  fx->SetParName(1, "Height1");
  fx->SetParName(2, "Centroid1");
  fx->SetParName(3, "Height2");
  fx->SetParName(4, "Centroid2");
  fx->SetParName(5, "Height3");
  fx->SetParName(6, "Centroid3");
  fx->SetParName(7, "Height4");
  fx->SetParName(8, "Centroid4");
  fx->SetParName(9, "BG offset");
  fx->SetParName(10,"BG slope");
  fx->SetParName(11,"BG Quad");

  fx->SetParameter(0, 5);
  fx->SetParameter(1, hist->GetBinContent(hist->FindBin(cent1)));
  fx->SetParameter(2, cent1);
  fx->SetParameter(3, hist->GetBinContent(hist->FindBin(cent2)));
  fx->SetParameter(4, cent2);
  fx->SetParameter(5, hist->GetBinContent(hist->FindBin(cent3)));
  fx->SetParameter(6, cent3);
  fx->SetParameter(7, hist->GetBinContent(hist->FindBin(cent4)));
  fx->SetParameter(8, cent4);
  fx->SetParameter(9, hist->GetBinContent(hist->FindBin(lower)));
  fx->SetParameter(10, -0.1);
  fx->SetParameter(11,bgquad);

  fx->SetParLimits(0,  0, 2);
  fx->SetParLimits(1,  0, (hist->GetMaximum()*10));
  fx->SetParLimits(3,  0, (hist->GetMaximum()*10));
  fx->SetParLimits(5,  0, (hist->GetMaximum()*10));
  fx->SetParLimits(7,  0, (hist->GetMaximum()*10));
  fx->SetParLimits(9,  -fabs(hist->GetMaximum()*10), fabs(hist->GetMaximum()*10));
  fx->SetParLimits(10, -100,100);
  fx->SetParLimits(11, -100,0);

  int npar = fx->GetNpar();
  if(bgquad == 0.0){
    fx->FixParameter((npar-1),0);
  }
  hist->Fit(fx,"","",lower,upper);
  
  double bgpar[3] = {fx->GetParameter(npar-3), fx->GetParameter(npar-2), fx->GetParameter(npar-1)};
  fbg->SetParameters(bgpar);
  fbg->SetLineColor(kBlack);
  fbg->SetLineStyle(9);
  f1->SetParameter(0, fx->GetParameter(0));
  f1->SetParameter(1, fx->GetParameter(1));
  f1->SetParameter(2, fx->GetParameter(2));
  f1->SetParameter(3, bgpar[0]);
  f1->SetParameter(4, bgpar[1]);
  f1->SetParameter(5, bgpar[2]);
  f1->SetLineColor(kBlue);
  f2->SetParameter(0, fx->GetParameter(0));
  f2->SetParameter(1, fx->GetParameter(3));
  f2->SetParameter(2, fx->GetParameter(4));
  f2->SetParameter(3, bgpar[0]);
  f2->SetParameter(4, bgpar[1]);
  f2->SetParameter(5, bgpar[2]);
  f2->SetLineColor(kGreen);
  f3->SetParameter(0, fx->GetParameter(0));
  f3->SetParameter(1, fx->GetParameter(5));
  f3->SetParameter(2, fx->GetParameter(6));
  f3->SetParameter(3, bgpar[0]);
  f3->SetParameter(4, bgpar[1]);
  f3->SetParameter(5, bgpar[2]);
  f3->SetLineColor(kMagenta);
  f4->SetParameter(0, fx->GetParameter(0));
  f4->SetParameter(1, fx->GetParameter(7));
  f4->SetParameter(2, fx->GetParameter(8));
  f4->SetParameter(3, bgpar[0]);
  f4->SetParameter(4, bgpar[1]);
  f4->SetParameter(5, bgpar[2]);
  f4->SetLineColor(kCyan);

  hist->GetListOfFunctions()->Add(f1);
  hist->GetListOfFunctions()->Add(f2);
  hist->GetListOfFunctions()->Add(f3);
  hist->GetListOfFunctions()->Add(f4);
  hist->GetListOfFunctions()->Add(fbg);

  double rebin = hist->GetBinWidth(hist->FindBin(lower));
  double bg    = fbg->Integral(lower,upper)/rebin;
  double sum   = hist->Integral(hist->FindBin(lower), hist->FindBin(upper)-1)-bg;
  double area  = fx->Integral(lower,upper)/rebin - bg;
  double area1 = f1->Integral(lower,upper)/rebin - bg;
  double area2 = f2->Integral(lower,upper)/rebin - bg;
  double area3 = f3->Integral(lower,upper)/rebin - bg;
  double area4 = f4->Integral(lower,upper)/rebin - bg;
  double chi2  = fx->GetChisquare()/(fx->GetNDF()-1);
	hist->GetListOfFunctions()->Add(fbg);
  system("Color 02");
  std::cout<< "Name: "  << "gausbg" << std::endl;
  std::cout<< "Area: "  << area  << std::endl;
  std::cout<< "Area1: " << area1 << std::endl;
  std::cout<< "Area2: " << area2 << std::endl;
  std::cout<< "Area3: " << area3 << std::endl;
  std::cout<< "Area4: " << area4 << std::endl;
  std::cout<< "Sum:  "  << sum  << std::endl;
  std::cout<< "Chi2: "  << chi2 << std::endl;
  gPad->Modified();
  gPad->Update();
}














