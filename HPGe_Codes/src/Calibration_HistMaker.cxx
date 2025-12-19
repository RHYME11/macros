//g++ Calibration_HistMaker.cxx -Wl,--no-as-needed `root-config --cflags --libs --glibs` -lSpectrum -lMinuit -lGuiHtml -lTreePlayer -lTMVA -L/opt/local/lib -lX11 -lXpm -O2 -Wl,--copy-dt-needed-entries -L/opt/local/lib -lX11 -lXpm `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -o Calibration_HistMaker


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
#include <TF1.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TSpectrum.h>
#include <Math/SpecFuncMathCore.h>
#include "TChannel.h"
#include "TTigress.h"
#include "TTigressHit.h"

TList *hlist;
std::vector<std::vector<double>> centroids(64);
std::vector<std::vector<double>> centroids_err(64);

// ================================ After this, need GRSISort Structure ======================== //
void Initialize(){
  hlist = new TList;
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

// ============================ Read Source File ========================================//
// Return vector<double> energies = energy of this source for the calibration
std::vector<double> ReadSourceFile(const std::string& source="60co"){
  std::vector<double> energies;
  std::vector<double> intensities;
  std::ifstream infile(Form("sources/%s.dat", source.c_str()));
  if (!infile.is_open()) {
    std::cerr << "Failed to open file: " << source << ".dat" << std::endl;
    return energies;
  }  
  std::string line;
  while (std::getline(infile, line)) {
    // Skip comment lines starting with '#'
    if (line.empty() || line[0] == '#') continue;
    std::stringstream ss(line);
    std::string isotope;
    double energy, intensity;
    // Only read the second column (after the isotope string)
    if (ss >> isotope >> energy >> intensity) {
      energies.push_back(energy);
      intensities.push_back(intensity);
    }
  }
  return energies;
}

// ============================ Read co60_linfit.dat File ========================================//
// Return gain and offset saved in the file. 
std::pair<std::vector<double>, std::vector<double>> ReadLinFitFile(const std::string& filename = "co60/co60_linfit.dat"){
  std::vector<double> lingain(64,-1.0);
  std::vector<double> linoffset(64,-1.0);
  std::ifstream infile(filename);
  if (!infile.is_open()) {
    std::cerr << "Failed to open file: " << filename << ".dat" << std::endl;
    return {lingain,linoffset};
  }  
  std::string line;
  while (std::getline(infile, line)) {
    // Skip comment lines starting with '#'
    if (line.empty() || line[0] == '#') continue;
    std::stringstream ss(line);
    int arryn;
    double g, o, fwhm, res;
    // Only read the second and third column (after the array number)
    if (ss >> arryn >> g >> o >> fwhm >> res) {
      lingain[arryn] = g;
      linoffset[arryn] = o;
    } // if over
  } // while loop over
  return {lingain, linoffset};
}

// ============================ Make the unclibrated energy ========================================//
// Make uncalibrated histogram
// Analysis TTree
void MakeRawHist(TChain *chain, char const *calfile){
  long nentries = chain->GetEntries();
  TTigress *tig = NULL;
  if(chain->FindBranch("TTigress")){
    chain->SetBranchAddress("TTigress", &tig);
  }else{
    std::cout << "Branch 'TTigress' not found! TTigress variable is NULL pointer" << std::endl;
    return;
  }
  
  TH1D *hs[64];
  for(int i=0;i<64;i++){
    hs[i] = new TH1D(Form("hs%i",i),Form("uncalibrated energy histogram at array %i",i), 4000,0,4000);
  }
  
  if(TChannel::ReadCalFile(calfile) < 1) {
    std::cout << "No channels found in calibration file " << calfile << "!" << std::endl;
    return;     
  }
  //TChannel::ReadCalFile(calfile); 
  std::cout<<std::endl;
  
  long xentry = 0;
  for(xentry;xentry<nentries;xentry++){
    chain->GetEntry(xentry);
    for(int i=0;i<tig->GetMultiplicity();i++){
      TTigressHit* tig_hit = tig->GetTigressHit(i);
      int arryn = tig_hit->GetArrayNumber(); // it will return xtal number, FulVA and FulVB from the same xtal will return same arraynumber. Eg, TIG01BN00A=0 TIG01GN00B=1 TIG05BN00A=16
      double charge = tig_hit->GetCharge();
      hs[arryn]->Fill(charge);
    }// loop xtal hits
    if((xentry%10000)==0){
      printf("Making Hist on entry: %lu / %lu \r", xentry, nentries);
      fflush(stdout);
    } 
  } // entries loop over 
  for(int i=0;i<64;i++){
    hlist->Add(hs[i]);
  }
  printf("Making Raw Hist DONE!  Entry: %lu / %lu \n", xentry, nentries);
}

// ============================ TF1: simple gaus + linear bg ===================================//
Double_t gaus_eqn(Double_t *x, Double_t *par){
  // Gussian part
  double A = par[0]; 
  double x0 = par[1];
  double sigma = par[2];
  double gaus = A*TMath::Exp(-(x[0]-x0)*(x[0]-x0))/(2*sigma*sigma);
  // linear bg
  double bg = par[3] + par[4]*x[0];
  
  return gaus + bg;
}

// ============================ Gaus Fit Raw Hist ========================================//
// this input type defined in the main function Step 5
void FitRawHist(const std::vector<std::vector<double>>& uncal_centroids){
  for(int i=0;i<64;i++){
    TH1D *hs = (TH1D *)hlist->FindObject(Form("hs%i",i));  
    if(hs->GetEntries()==0) continue;
    for(int j=0;j<uncal_centroids[i].size();j++){
      Int_t bin_guess = hs->FindBin(uncal_centroids[i][j]);
      hs->SetAxisRange(bin_guess-15, bin_guess+15, "X");
      Int_t peak_bin = hs->GetMaximumBin();
      double x_guess = hs->GetBinCenter(peak_bin);
      double y_guess = hs->GetBinContent(peak_bin);
      hs->GetXaxis()->SetRange(0,0); // unzoom
      TF1 *fx = new TF1(Form("fx%i_peak%i",i,j), gaus_eqn, x_guess-8, x_guess+8,5);
      fx->SetParameters(y_guess, x_guess, 0.5, y_guess/100., -0.1);
      fx->SetParLimits(0, y_guess*0.8, y_guess*1.2); //area
      fx->SetParLimits(1, x_guess - 10, x_guess + 10); //centroid
      fx->SetParLimits(2, 0., 15); //sigma of gaussian distribution
      //fx->SetParLimits(5, -10, -0.1); //background noise constant
      hs->Fit(fx, "QR+");
      centroids[i].push_back(fx->GetParameter(1));           
      centroids_err[i].push_back(fx->GetParError(1));           
    } // j (peak number) loop over
  } //i (array number) loop over
}

// ====================================== main() ==========================================//
// argv1: CalibrationFile
// argv2: Source Name
// argv3...: AnalysisTree File Path
int main(int argc, char** argv){

  if(argc<4){
    printf("Input Calibration file, source and Analysistree file paths\n");
    return 1;
  }
  //Step 1: loop over root file if files are valid
  TChain *chain = new TChain("AnalysisTree");
  for(int i=3;i<argc;i++){
    std::string rootfilename = argv[i];
    chain->Add(rootfilename.c_str());
  }
  if(chain->GetEntries()==0){
    printf("No valid root file input\n");
    return 1;
  }
  
  // Step 2: make uncalibrated energy histogram
  char const *calfile = argv[1];
  Initialize();
  MakeRawHist(chain, calfile);  
  
  // Step 3: Readout co60_linfit.dat gain and offset
  auto [lingain, linoff] = ReadLinFitFile();
  
  // Step 4: handle source, read out energes for calibration 
  std::string source = FormatIsotopeName(argv[2]);
  std::vector<double> energies = ReadSourceFile(source);
  
  // Step 5: calculated uncalibrated energy centroid of peaks based on Step 3 and Step 4
  // uncal_centroids[i][j] = centroid in the uncalibrated hs_i(array number) related #j energy from source.dat
  std::vector<std::vector<double>> uncal_centroids(64, std::vector<double>(energies.size()));
  for(int i=0;i<lingain.size();i++){
    if(lingain[i]<0) continue; // skip non-existent crystals
    for(int j=0;j<energies.size();j++){
      // uncal_E = (cal_E - offset) / gain
      uncal_centroids[i][j] = (energies[j]-linoff[i])/lingain[i];
      if(i==38) printf("j = %i, energies = %.2f, uncal = %.2f\n",j, energies[j], uncal_centroids[i][j]);
    }// j (ref energy) loop over
  } // i (array number) loop over

  // Step 6: Only fit the peak but don't fit the calibration. Return uncalibrated centroids(=centroids[])
  FitRawHist(uncal_centroids);
  
  // Step 7: Write array_number, uncal_e and energies to file, and save hlist into root file
  std::ofstream outfile(Form("peaks_%s.dat",source.c_str()));
  outfile << "# ArrayNum\tUncal_E\tErr\tEnergies\n";
  for(int i=0; i<64; i++){
    if(centroids[i].size() == 0) continue; // skip empty array
    for(int j=0; j<centroids[i].size(); j++){
      outfile << i << '\t'
              << centroids[i][j] << '\t'
              << centroids_err[i][j] << '\t'
              << energies[j] << '\n';
    }
  }
  outfile.close(); 
 
  TFile *newf = new TFile(Form("peaks_%s.root",source.c_str()), "recreate");
  newf->cd();
  hlist->Write();
  newf->Close(); 

  return 0;
}


















