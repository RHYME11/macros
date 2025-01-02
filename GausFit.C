
// ==== Single ==== // 
double Gaus(double *dim, double *par){
  // 3 parameters: a, x0, sigma;
  double a = par[0];
  double x0 = par[1];
  double sigma = par[2];
  double x = dim[0];

  return a*TMath::Exp(-((x-x0)*(x-x0))/(2*sigma*sigma));
}

double GausBG(double *dim, double *par){
  // 3 parameters: a, x0, sigma;
  double a = par[0];
  double x0 = par[1];
  double sigma = par[2];
  double x = dim[0];

  return a*TMath::Exp(-((x-x0)*(x-x0))/(2*sigma*sigma)) + (par[3]*x + par[4]);
}

// ==== Double ==== // 
double DoubleGaus(double *dim, double *par){
  double a1 = par[0];
  double x1 = par[1];
  double a2 = par[2];
  double x2 = par[3];
  double sigma = par[4];
  double x = dim[0];

  double arr1[3] = {a1,x1,sigma};
  double arr2[3] = {a2,x2,sigma};
  return Gaus(dim,arr1) + Gaus(dim,arr2);
}

double DoubleGausBG(double *dim, double *par){
  double a1 = par[0];
  double x1 = par[1];
  double a2 = par[2];
  double x2 = par[3];
  double sigma = par[4];
  double x = dim[0];

  double arr1[3] = {a1,x1,sigma};
  double arr2[3] = {a2,x2,sigma};
  return Gaus(dim,arr1) + Gaus(dim,arr2) + (par[5]*x + par[6]);
}

// ==== Triple ==== // 
double TripleGaus(double *dim, double *par){
  double a1 = par[0];
  double x1 = par[1];
  double a2 = par[2];
  double x2 = par[3];
  double a3 = par[4];
  double x3 = par[5];
  double sigma = par[6];
  double x = dim[0];

  double arr1[3] = {a1,x1,sigma};
  double arr2[3] = {a2,x2,sigma};
  double arr3[3] = {a3,x3,sigma};
  return Gaus(dim,arr1) + Gaus(dim,arr2) + Gaus(dim,arr3);
}

double TripleGausBG(double *dim, double *par){
  double a1 = par[0];
  double x1 = par[1];
  double a2 = par[2];
  double x2 = par[3];
  double a3 = par[4];
  double x3 = par[5];
  double sigma = par[6];
  double x = dim[0];

  double arr1[3] = {a1,x1,sigma};
  double arr2[3] = {a2,x2,sigma};
  double arr3[3] = {a3,x3,sigma};
  return Gaus(dim,arr1) + Gaus(dim,arr2) + Gaus(dim,arr3) + (par[7]*x + par[8]);
}



// ==== Fitting Functions ==== //
void GausFit(TH1 *hist, double lower, double upper){
  hist->GetListOfFunctions()->Clear();

  TF1 *fx = new TF1("fx", GausBG, lower, upper, 5);
  TF1 *fbg= new TF1("fbg","[0]*x+[1]",lower, upper);

  fx->SetParName(0, "Height");
  fx->SetParName(1, "Centroid");
  fx->SetParName(2, "Sigma");
  fx->SetParName(3, "Bg slope");
  fx->SetParName(4, "Bg offset");

  fx->SetParameter(0, hist->GetMaximum());
  fx->SetParameter(1, (lower+upper)/2.);
  fx->SetParameter(2, 5);
  fx->SetParameter(3, -0.1);
  fx->SetParameter(4, hist->GetMinimum());

  fx->SetParLimits(0,0,1e6);
  fx->SetParLimits(3,-100,1.);
  fx->SetParLimits(4,-1e6,1e6);

  hist->Fit(fx,"","",lower,upper);
  fbg->SetParameter(0, fx->GetParameter(3));
  fbg->SetParameter(1, fx->GetParameter(4));
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

void DoubleGausFit(TH1 *hist, double cent1, double cent2,double lower, double upper){
  hist->GetListOfFunctions()->Clear();

  TF1 *fx = new TF1("fx", DoubleGausBG, lower, upper, 7);
  TF1 *f1 = new TF1("f1", GausBG, lower, upper, 5);
  TF1 *f2 = new TF1("f2", GausBG, lower, upper, 5);
  TF1 *fbg= new TF1("fbg","[0]*x+[1]",lower, upper);

  fx->SetParName(0, "Height 1");
  fx->SetParName(1, "Centroid 1");
  fx->SetParName(2, "Height 2");
  fx->SetParName(3, "Centroid 2");
  fx->SetParName(4, "Sigma");
  fx->SetParName(5, "Bg slope");
  fx->SetParName(6, "Bg offset");

  fx->SetParameter(0, hist->GetMaximum());
  fx->SetParameter(1, cent1);
  fx->SetParameter(2, hist->GetMaximum());
  fx->SetParameter(3, cent2);
  fx->SetParameter(4, 1);
  fx->SetParameter(5, -0.1);
  fx->SetParameter(6, hist->GetMinimum());

  fx->SetParLimits(0,0,.1e6);
  fx->SetParLimits(2,0,.1e6);
  fx->SetParLimits(5,-100,.1);
  fx->SetParLimits(6,-1e6,1e6);
	
	hist->Fit(fx,"","",lower,upper);
  f1->SetParameter(0, fx->GetParameter(0));
  f1->SetParameter(1, fx->GetParameter(1));
  f1->SetParameter(2, fx->GetParameter(4));
  f1->SetParameter(3, fx->GetParameter(5));
  f1->SetParameter(4, fx->GetParameter(6));
  f1->SetLineColor(kBlue);
  f2->SetParameter(0, fx->GetParameter(2));
  f2->SetParameter(1, fx->GetParameter(3));
  f2->SetParameter(2, fx->GetParameter(4));
  f2->SetParameter(3, fx->GetParameter(5));
  f2->SetParameter(4, fx->GetParameter(6));
  f2->SetLineColor(kGreen);
  fbg->SetParameter(0, fx->GetParameter(5));
  fbg->SetParameter(1, fx->GetParameter(6));
  fbg->SetLineColor(kBlack);
  fbg->SetLineStyle(9);

  double rebin = hist->GetBinWidth(hist->FindBin(lower));
  double bg = fbg->Integral(lower,upper)/rebin;
  double sum  = hist->Integral(hist->FindBin(lower), hist->FindBin(upper)-1)-bg;
  double area  = fx->Integral(lower,upper)/rebin - bg;
  double area1 = f1->Integral(lower,upper)/rebin - bg;
  double area2 = f2->Integral(lower,upper)/rebin - bg;
  double chi2 = fx->GetChisquare()/(fx->GetNDF()-1);

  hist->GetListOfFunctions()->Add(f1);
  hist->GetListOfFunctions()->Add(f2);
  hist->GetListOfFunctions()->Add(fbg);
  system("Color 02");
  std::cout<< "Name: " << "gausbg" << std::endl;
  std::cout<< "Area: " << area  << std::endl;
  std::cout<< "Area1: " << area1  << std::endl;
  std::cout<< "Area2: " << area2  << std::endl;
  std::cout<< "Sum:  " << sum  << std::endl;
  std::cout<< "Chi2: " << chi2 << std::endl;
  gPad->Modified();
  gPad->Update();
}


void TripleGausFit(TH1 *hist, double cent1, double cent2, double cent3, double lower, double upper){
  hist->GetListOfFunctions()->Clear();

  TF1 *fx = new TF1("fx", TripleGausBG, lower, upper, 9);
  TF1 *f1 = new TF1("f1", GausBG, lower, upper, 5);
  TF1 *f2 = new TF1("f2", GausBG, lower, upper, 5);
  TF1 *f3 = new TF1("f3", GausBG, lower, upper, 5);
  TF1 *fbg= new TF1("fbg","[0]*x+[1]",lower, upper);

  fx->SetParName(0, "Height 1");
  fx->SetParName(1, "Centroid 1");
  fx->SetParName(2, "Height 2");
  fx->SetParName(3, "Centroid 2");
  fx->SetParName(4, "Height 3");
  fx->SetParName(5, "Centroid 3");
  fx->SetParName(6, "Sigma");
  fx->SetParName(7, "Bg slope");
  fx->SetParName(8, "Bg offset");

  fx->SetParameter(0, hist->GetMaximum());
  fx->SetParameter(1, cent1);
  fx->SetParameter(2, hist->GetMaximum());
  fx->SetParameter(3, cent2);
  fx->SetParameter(4, hist->GetMaximum());
  fx->SetParameter(5, cent3);
  fx->SetParameter(6, 1);
  fx->SetParameter(7, -0.1);
  fx->SetParameter(8, hist->GetMinimum());

  fx->SetParLimits(0,0,(hist->GetMaximum())*1.5);
  fx->SetParLimits(2,0,(hist->GetMaximum())*1.5);
  fx->SetParLimits(4,0,(hist->GetMaximum())*1.5);
  fx->SetParLimits(7,-100,.1);
  fx->SetParLimits(8,-1e6,1e6);
	
	hist->Fit(fx,"","",lower,upper);
  f1->SetParameter(0, fx->GetParameter(0));
  f1->SetParameter(1, fx->GetParameter(1));
  f1->SetParameter(2, fx->GetParameter(6));
  f1->SetParameter(3, fx->GetParameter(7));
  f1->SetParameter(4, fx->GetParameter(8));
  f1->SetLineColor(kBlue);
  f2->SetParameter(0, fx->GetParameter(2));
  f2->SetParameter(1, fx->GetParameter(3));
  f2->SetParameter(2, fx->GetParameter(6));
  f2->SetParameter(3, fx->GetParameter(7));
  f2->SetParameter(4, fx->GetParameter(8));
  f2->SetLineColor(kGreen);
  f3->SetParameter(0, fx->GetParameter(4));
  f3->SetParameter(1, fx->GetParameter(5));
  f3->SetParameter(2, fx->GetParameter(6));
  f3->SetParameter(3, fx->GetParameter(7));
  f3->SetParameter(4, fx->GetParameter(8));
  f3->SetLineColor(kMagenta);
  fbg->SetParameter(0, fx->GetParameter(7));
  fbg->SetParameter(1, fx->GetParameter(8));
  fbg->SetLineColor(kBlack);
  fbg->SetLineStyle(9);

  double rebin = hist->GetBinWidth(hist->FindBin(lower));
  double bg = fbg->Integral(lower,upper)/rebin;
  double sum  = hist->Integral(hist->FindBin(lower), hist->FindBin(upper)-1)-bg;
  double area  = fx->Integral(lower,upper)/rebin - bg;
  double area1 = f1->Integral(lower,upper)/rebin - bg;
  double area2 = f2->Integral(lower,upper)/rebin - bg;
  double area3 = f3->Integral(lower,upper)/rebin - bg;
  double chi2 = fx->GetChisquare()/(fx->GetNDF()-1);

  hist->GetListOfFunctions()->Add(f1);
  hist->GetListOfFunctions()->Add(f2);
  hist->GetListOfFunctions()->Add(f3);
  hist->GetListOfFunctions()->Add(fbg);
  system("Color 02");
  std::cout<< "Name: " << "gausbg" << std::endl;
  std::cout<< "Area: " << area  << std::endl;
  std::cout<< "Area1: " << area1  << std::endl;
  std::cout<< "Area2: " << area2  << std::endl;
  std::cout<< "Area3: " << area3  << std::endl;
  std::cout<< "Sum:  " << sum  << std::endl;
  std::cout<< "Chi2: " << chi2 << std::endl;
  gPad->Modified();
  gPad->Update();
}



