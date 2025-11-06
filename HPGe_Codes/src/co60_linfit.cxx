//g++ co60_linfit.cxx -Wl,--no-as-needed `root-config --cflags --libs --glibs` -lSpectrum -lMinuit -lGuiHtml -lTreePlayer -lTMVA -L/opt/local/lib -lX11 -lXpm -O2 -Wl,--copy-dt-needed-entries -L/opt/local/lib -lX11 -lXpm `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -o co60_linfit


#include <iostream>
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
TList *glist;
double gain[64];
double offset[64];
double sigmas[64];

// ================================ After this, need GRSISort Structure ======================== //
void Initialize(){
  hlist = new TList;
  glist = new TList; 
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


// =========================== Write information into file ================================ //
void WriteFile(){
  std::ofstream outfile("co60_linfit.dat");
  outfile << "# ArrayNum  Gain   Offset  FWHM(1332)  Res%(1332)\n";
  for (int i = 0; i < 64; ++i) {
    if (gain[i] == 0 && offset[i] == 0 && sigmas[i] == 0) continue; // skip unused channels
    double fwhm = sigmas[i] * 2.35;
    double resolution = fwhm / 1332.0 * 100;
    outfile << std::setw(9) << i
            << std::setw(10) << std::fixed << std::setprecision(6) << gain[i]
            << std::setw(10) << offset[i]
            << std::setw(13) << fwhm
            << std::setw(11) << resolution << '\n';
  }
  outfile.close();
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
    hs[i] = new TH1D(Form("hs%i",i),Form("unclibrated energy histogram at array %i",i), 4000,0,4000);
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

// ============================ PeakHunt ====================================//
// Peak search in TH1 *hist by using TSpectrum;
// Return a vector of x-axis position of nref peaks with highest y-values;
// TH1 *hist needs to zoom in a reasonable range to skip noise at low ADC channel. 
std::vector<Double_t> PeakHunt(TH1 *hist, int nref=2){
  TSpectrum *s = new TSpectrum(10); //max positions = 10
  hist->GetXaxis()->SetRangeUser(20,1e4); // skip the nosiy peak
  Int_t npeaks = s->Search(hist,2,"",0.13); // rough peak-search through the entire hist
  Double_t *xpeaks = s->GetPositionX();
  Double_t *ypeaks = s->GetPositionY();
  std::vector<std::pair<Double_t, Double_t>> peaks;
  for(int ipeak=0;ipeak<npeaks;ipeak++){
    peaks.emplace_back(xpeaks[ipeak], ypeaks[ipeak]);
  }  
  // order peaks based on y-values from highest to lowest;
  // extract the first nref elements (nref highest peak) and put their x-values in a new vector top_xpeaks;
  // then order top_xpeaks from lowest to highest based on their values;
  // extract the first and last bin, which should be centers of the first and last true alpha peaks;
  std::sort(peaks.begin(), peaks.end(), [](const auto &m, const auto &n){return m.second > n.second;});
  std::vector<Double_t> top_xpeaks;
  if(npeaks>=nref){
    for(int ipeak=0;ipeak<nref;ipeak++){
      top_xpeaks.push_back(peaks[ipeak].first);
    }
  }else{
    for(int ipeak=0;ipeak<npeaks;ipeak++){
      top_xpeaks.push_back(peaks[ipeak].first);
    }
  }  
  std::sort(top_xpeaks.begin(), top_xpeaks.end());
  return top_xpeaks;
}

// ============================ TF1: fancy gaus fitting function ========================================//
// I don't understand "erf" part;
// I copied from Stephen's code;
// He probably copied from Red-Ware;
//fit[j] = new TF1(Form("Fit %i-%i", i, j), "[0]*(exp(-((x-[1])^2/(2*[2]^2))))*(1+ROOT::Math::erf([5]*((x-[1]))/([2]*pow(2,0.5)))) + [3] + [4]*(0.5*(1-(ROOT::Math::erf(((x-[1])/([2]*2^(0.5)))))))", x_pos - 50, x_pos + 50);
//fit[j]->SetParameters(y_pos, x_pos, 1, 15, 1, -1);
//fit[j]->SetParLimits(0, 10, 1e6); //area
//fit[j]->SetParLimits(1, x_pos - 10, x_pos + 10); //centroid
//fit[j]->SetParLimits(2, 0.2, 15); //sigma
//fit[j]->SetParLimits(4, 0.1, 100); //magnitude of step in background noise
//fit[j]->SetParLimits(5, -10, -0.1); //background noise constant
Double_t peak_eqn(Double_t *x, Double_t *par){
 
  Double_t gaus_val = par[0]*exp(-(pow((x[0]-par[1]),2))/(2*par[2]*par[2]));
  Double_t erf_inp  = (x[0]-par[1])/(par[2]*pow(2,0.5));
  Double_t return_val = gaus_val * (1+ROOT::Math::erf(par[5]*erf_inp)) 
                      + par[3] + par[4]*(0.5*(1-ROOT::Math::erf(erf_inp)));
  return return_val;
}

// ============================ CalRawHist(): Fit 60Co ====================================//
void CalRawHist(std::vector<double> energies, std::vector<double> energy_err){
  int nref = energies.size(); // how many peaks used for the calibration (nref = 2 for 60Co)
  for(int i=0;i<64;i++){
    TH1D *hs = (TH1D *)hlist->FindObject(Form("hs%i",i));
    if(hs->GetEntries()==0) continue;
    std::vector<Double_t> xpeaks = PeakHunt(hs, nref);
    if(xpeaks.size()<nref){
      printf("Arraynumber[%i] has %i peaks less than %i peaks listed in source.dat", i, xpeaks.size(), nref); 
    }else{
      double centroids[2];
      double centroid_errs[2];
      double sigma = -1;
      for(int j=0;j<xpeaks.size();j++){
        double binc = hs->GetBinContent(hs->FindBin(xpeaks[j]));
        TF1 *fx = new TF1(Form("fx%i_peak%i",i,j), peak_eqn, xpeaks[j]-50, xpeaks[j]+50,5);
        fx->SetParameters(binc, xpeaks[j], 1, 15, 1, -1);
        fx->SetParLimits(0, 10, 1e6); //area
        fx->SetParLimits(1, xpeaks[j] - 10, xpeaks[j] + 10); //centroid
        fx->SetParLimits(2, 0.2, 15); //sigma
        fx->SetParLimits(4, 0.1, 100); //magnitude of step in background noise
        fx->SetParLimits(5, -10, -0.1); //background noise constant
        hs->Fit(fx,"RQ+");
        centroids[j] = fx->GetParameter(1);
        centroid_errs[j] = fx->GetParError(1);
        sigma = fx->GetParameter(2); // Report resolution with sigma from 1332-keV peak
      }// peaks loop over
      TGraphErrors *gr = new TGraphErrors(nref, centroids, energies.data(), centroid_errs, energy_err.data());
      gr->SetName(Form("gr%i",i));
      glist->Add(gr);
      TF1 *flin = new TF1(Form("flin%i",i),"[0]+[1]*x");
      flin->SetParameters(10,1);
      gr->Fit(flin, "Q+"); 
      offset[i] = flin->GetParameter(0);
      gain[i] = flin->GetParameter(1);
      sigmas[i] = sigma * gain[i];
    }// else 
  } // hist loop over   
}



// ====================================== main() ==========================================//
// argv1: CalibrationFile
// argv2...: AnalysisTree File Path
int main(int argc, char** argv){

  if(argc<3){
    printf("Input Calibration file and Analysistree file path");
    return 1;
  }
  //Step 1: loop over root file if files are valid
  TChain *chain = new TChain("AnalysisTree");
  for(int i=2;i<argc;i++){
    std::string rootfilename = argv[i];
    chain->Add(rootfilename.c_str());
  }
  if(chain->GetEntries()==0){
    printf("No valid root file input\n");
    return 1;
  }
  // Step 2: make uncalibrated histogram
  char const *calfile = argv[1];
  Initialize();
  MakeRawHist(chain, calfile);  
   

  // Step 3: handle source
  std::string source = FormatIsotopeName(); // default 60co
  std::vector<double> energies = ReadSourceFile(source);
  std::vector<double> energy_err(energies.size(), 0.0);  // when calibrate, no error on gam energy
  
  // Step 4: Calibrate 60Co
  CalRawHist(energies, energy_err); 
  
  // Step 5: Write gain, offset and sigma into file, and save hlist into root file
  WriteFile(); 
  TFile *newf = new TFile("co60_linfit.root", "recreate");
  newf->cd();
  hlist->Write();
  glist->Write();
  newf->Close(); 

  return 0;
}


















