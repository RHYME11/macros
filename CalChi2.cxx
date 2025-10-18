//g++ CalChi2.cxx -o bin/CalChi2 `root-config --cflags --glibs`

#include <iostream>
#include <vector>
#include <regex>
#include <string>
#include <cmath>
#include <iomanip>


#include <TFile.h>
#include <TH1.h>
#include <TGraph.h>
#include <TCanvas.h>

//========================================================================//
Double_t Chi2(double obs, double exp, double sigma = -1){
  if(sigma<0){
    sigma = sqrt(obs);
  }
  return pow((exp-obs),2)/pow(sigma,2);
}

//========================================================================//
Double_t Chi2_LH(double obs, double exp){
  return 2*(exp - obs + obs*log(obs/exp));
}

//========================================================================//
Double_t CalChi2(TH1D *hdat, TH1D *hsim, double lowerE, double upperE, double scal=1.){
  int lowerbin =  hdat->FindBin(lowerE);
  int upperbin =  hdat->FindBin(upperE);
  double chi2 = 0;
  for(int binnum=lowerbin;binnum<upperbin;binnum++){
    double obs = hdat->GetBinContent(binnum);
    double exp = hsim->GetBinContent(binnum);
    if(obs==0 || exp==0) continue;
    exp = exp * scal;
    //chi2 += Chi2_LH(obs, exp);
    chi2 += Chi2(obs, exp);
  }
  return chi2;
}





//=== main function ===//
// Input1: lower E (not bin number);
// Input2: upper E (not bin number);
// Input3: data root file
// Input4: simulation root file 
int main(int argc, char **argv){

  if(argc<2){
    printf("Add Inputs! \n");
    return 1;
  }
  if(argc<4){
    printf("Add ROOT files!\n");
    return 1;
  }

  TFile *infile_dat = TFile::Open(argv[3]);
  TFile *infile_sim = TFile::Open(argv[4]);
  
  TH1D *hdat = (TH1D *)infile_dat->Get("hac6");
  TH1D *hsim = (TH1D *)infile_sim->Get("hac6");
  hdat->SetName("hdat");
  hsim->SetName("hsim");  
 
  hdat->Rebin(10);
  hsim->Rebin(10); 
 
  double lowerE = std::stod(argv[1]); 
  double upperE = std::stod(argv[2]); 
  
  for(double fac=0.8; fac<1.2; fac=fac*1.001){
    std::cout << std::fixed << std::setprecision(4)
              << fac        << "\t"
              << CalChi2(hdat, hsim, lowerE, upperE, fac) << std::endl;
  }

  return 0;
}




