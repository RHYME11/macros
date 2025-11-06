//g++ Calibration.cxx -Wl,--no-as-needed `root-config --cflags --libs --glibs` -lSpectrum -lMinuit -lGuiHtml -lTreePlayer -lTMVA -L/opt/local/lib -lX11 -lXpm -O2 -Wl,--copy-dt-needed-entries -L/opt/local/lib -lX11 -lXpm `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -o Calibration


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cctype>
#include <algorithm>
#include <iomanip>
#include <cmath>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TSpectrum.h>
#include <Math/SpecFuncMathCore.h>
#include "TChannel.h"
#include "TTigress.h"
#include "TTigressHit.h"

TList *glist;
TH2D *sumc;
TH2D *sume;
double quad[64];
double gain[64];
double offset[64];
std::vector<std::vector<double>> uncalE(64);
std::vector<std::vector<double>> uncalE_err(64);
std::vector<std::vector<double>> energies(64);
std::vector<std::vector<double>> energies_err(64);
// ================================ After this, need GRSISort Structure ======================== //
void Initialize(){
  glist = new TList; 
  sumc = new TH2D("sumc","uncalibrated summary plot", 4000,0,4000,64,0,64);
  sume = new TH2D("sume","calibrated summary plot"  , 4000,0,4000,64,0,64);
}
// ============================ Source Name Format ========================================//
// convert input source name to number + lowercase letter (eg, 60co, 133ba)
std::string FormatIsotopeName(const std::string& srcName="60co") {
  std::string number, letters;
  for (char c : srcName) {
    if (std::isdigit(c))
      number += c;
    else if (std::isalpha(c))
      letters += std::tolower(c);  // ensure lowercase                                                                                                                                                             
  }       
  return number + letters;
}

// ====================================== ReadPeaksFile() ========================================//
void ReadPeaksFile(const std::string& source){
  std::ifstream infile(Form("peaks/peaks_%s.dat", source.c_str()));
  if (!infile.is_open()) {
    std::cerr << "Cannot open dat file: peaks_" << source << ".dat" << std::endl;
    return;
  }
  std::string line;
  while(std::getline(infile, line)){
    if(line.empty() || line[0] == '#') continue;    
    std::stringstream ss(line);
    int arryn;
    double uncal_E, err, energy;
    if(ss >> arryn >> uncal_E >> err >> energy){
      uncalE[arryn].push_back(uncal_E);
      uncalE_err[arryn].push_back(err);
      energies[arryn].push_back(energy);
    }
  }    
}

// ====================================== Calibrate() ========================================//
void WriteCalibFile(const std::string& filename="cal_pars.dat") {
  std::ofstream fout(filename);
  fout << "float non_lin[64] = {";
  for (int i = 0; i < 64; ++i) {
    fout << quad[i];
    if (i != 63) fout << ", ";
  }
  fout << "};\n";

  fout << "float gain[64] = {";
  for (int i = 0; i < 64; ++i) {
    fout << gain[i];
    if (i != 63) fout << ", ";
  }
  fout << "};\n";

  fout << "float offset[64] = {";
  for (int i = 0; i < 64; ++i) {
    fout << offset[i];
    if (i != 63) fout << ", ";
  }
  fout << "};\n";

  fout.close();
}

// ====================================== Calibrate() ========================================//
void Calibrate(){
  for(int i=0;i<uncalE.size();i++){
    std::cout << "Fitting array " << i 
              << ": N points = " << uncalE[i].size() << std::endl;
    if(uncalE[i].empty()) continue;
    energies_err[i].resize(energies[i].size(), 0.0);
    TGraphErrors *gr = new TGraphErrors(uncalE[i].size(),
                                        uncalE[i].data(), energies[i].data(), 
                                        uncalE_err[i].data(),energies_err[i].data());
    gr->SetName(Form("gr%i",i));
    TF1 *fx = new TF1(Form("fx%i",i), "[0]+[1]*x+[2]*x*x");
    fx->SetParameters(30,1.5,1e-6);
    gr->Fit(fx,"Q");
    offset[i] = fx->GetParameter(0);
    gain[i] = fx->GetParameter(1);
    quad[i] = fx->GetParameter(2);
    glist->Add(gr);
  } // i (array number) loop over
}

// ====================================== DrawSum(): Draw Summary Plot ====================================//
void DrawSum(std::vector<std::string> sources){
  // arraynum -> (bin center) -> content
  std::map<int, std::map<double, double>> unc; 
  std::map<int, std::map<double, double>> cal; 
  for(int i=0;i<sources.size();i++){
    TFile *curfile = TFile::Open(Form("peaks_%s.root",sources[i].c_str()));
    if (!curfile || curfile->IsZombie()) {
      std::cerr << "Cannot open file: " << sources[i] << std::endl;
      continue;
    }
    for(int j=0;j<64;j++){
      TH1D *htemp = (TH1D *)curfile->Get(Form("hs%i",j));
      if(! htemp || htemp->GetEntries()==0) continue;
      for(int n=1;n<=htemp->GetNbinsX();n++){
        double charge = htemp->GetBinCenter(n);
        double energy = offset[j] + gain[j]*charge + quad[j]*charge*charge;
        double content = htemp->GetBinContent(n);
        unc[j][charge] += content;
        cal[j][energy] += content;
      }// n (# of bin in the histogram) loop over
    } // j (array number) loop over
    curfile->Close();
  } // i (sources) loop over 
  
  for (const auto& [arraynum, charge_map] : unc) {
    for (const auto& [charge, content] : charge_map) {
      sumc->Fill(charge, arraynum, content);  
    }
  }
  for (const auto& [arraynum, energy_map] : cal) {
    for (const auto& [energy, content] : energy_map) {
      sume->Fill(energy, arraynum, content);  
    }
  }
}

// ====================================== main() ==========================================//
// argv1...: sources name
int main(int argc, char** argv){

  Initialize();
  
  if(argc<2){
    printf("Input source name");
    return 1;
  }
  // Step 1: handle input argv
  std::vector<std::string> sources;
  for(int i=1;i<argc;i++){
    sources.push_back(FormatIsotopeName(argv[i])); 
  }

  // Step 2: readout all peaks_source.dat and put them into 2D vectors
  for(int i=0;i<sources.size();i++){
    ReadPeaksFile(sources[i]);  
  }
  
  // Step 3: Calibration 
  Calibrate();

  // Step 4: Draw summary plot
  DrawSum(sources); 
  
  // Step 5: write calibration coefficiency, and write the root file
  WriteCalibFile();
  TFile *newf = new TFile("calibration.root", "recreate");
  newf->cd();
  sumc->Write();
  sume->Write();
  glist->Write();
  newf->Close(); 

  return 0;
}


















