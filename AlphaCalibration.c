//g++ AlphaCalibration.c -Wl,--no-as-needed `root-config --cflags --libs --glibs` -lSpectrum -lMinuit -lGuiHtml -lTreePlayer -lTMVA -L/opt/local/lib -lX11 -lXpm -O2 -Wl,--copy-dt-needed-entries -L/opt/local/lib -lX11 -lXpm `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -o bin/Alphacal






#include <iostream>  
#include <iomanip>   
#include <fstream>   
#include <cmath>     
#include <algorithm>
#include <map>
#include <utility>
#include <string> 
#include <stdio.h>   
#include <stdlib.h>  
#include "TCutG.h"   
#include "TH1.h"     
#include "TF1.h"     
#include "TTree.h"   
#include "TChain.h"  
#include "TH2.h"     
#include "TFile.h"   
#include "TDirectory.h"
#include "TList.h"   
#include "TTip.h"    
#include "TTigress.h"
#include "TRF.h"     
#include "TSpectrum.h"
#include "TChannel.h"
#include "TParserLibrary.h"
#include "TEnv.h"    
#include <dirent.h>

TList *hlist;
TList *flist;
TH2D *sumc;
TH2D *sume;
std::map<int,std::pair<double, double>> calmap; // <channel, <gain,offset>>


// ============= TripleAlphaHighE_Fun() ======================= //
// TF1 formula for triple alpha source, 239Pu, 241Am, and 244Cm
// Simulate a triple alpha spectrum for comparison to a histogram
// Energy calibration and width of peaks are parameters given up to quadratic linear
Double_t TripleAlphaHighE_Fun(Double_t *x, Double_t *par){
// Parameters:
// par[0] -- Normalization factor for 239Pu group
// par[1] -- Normalization factor for 241Am group
// par[2] -- Normalization factor for 244Cm group
// par[3] -- FWHM of peaks in keV
// par[4] -- Bg: Constant offset in counts*
// par[5] -- Width of Pu peaks relative to Cm
// par[6] -- Width of Am peaks relative to Cm
// par[7] -- Offset
// par[8] -- Linear gain
  
  Double_t E = par[7]+par[8]*x[0];
  Double_t sigmaCm = par[3]/2.35;
  Double_t sigmaPu = sigmaCm*par[5];
  Double_t sigmaAm = sigmaCm*par[6];

  Double_t return_f = par[0] * ( 0.7077 * TMath::Gaus(E,5156.59,sigmaPu)    // Pu
                              +0.1711 * TMath::Gaus(E,5144.30,sigmaPu)
                              +0.1194 * TMath::Gaus(E,5105.80,sigmaPu) )
                     +par[1] * ( 0.0036 * TMath::Gaus(E,5544.5,sigmaAm)     // Am
                              +0.0166 * TMath::Gaus(E,5388,sigmaAm)
                              +0.848 * TMath::Gaus(E,5485.56,sigmaAm)
                              +0.131 * TMath::Gaus(E,5442.8,sigmaAm) )
                     +par[2] * ( 0.769 * TMath::Gaus(E,5804.77,sigmaCm)     // Cu
                              +0.231 * TMath::Gaus(E,5762.16,sigmaCm ) )
                     +par[4];                                               // Bg
  return return_f;
}


// ============= TripleAlphaLowE_Fun() ======================= //
// TF1 formula for triple alpha source, 148Gd, 230Th, 244Cm
// Simulate a triple alpha spectrum for comparison to a histogram
// Energy calibration and width of peaks are parameters given up to quadratic linear
Double_t TripleAlphaLowE_Fun(Double_t *x, Double_t *par){
// Parameters:
// par[0] -- Normalization factor for 148Gd group
// par[1] -- Normalization factor for 230Th group
// par[2] -- Normalization factor for 244Cm group
// par[3] -- FWHM of peaks in keV
// par[4] -- Bg: Constant offset in counts*
// par[5] -- Width of Gd peaks relative to Cm
// par[6] -- Width of Th peaks relative to Cm
// par[7] -- Offset
// par[8] -- Linear gain
  
  Double_t E = par[7]+par[8]*x[0];
  Double_t sigmaCm = par[3]/2.35;
  Double_t sigmaGd = sigmaCm*par[5];
  Double_t sigmaTh = sigmaCm*par[6];

  Double_t return_f = par[0] * ( 1.0 * TMath::Gaus(E,3182.690,sigmaGd))
                     +par[1] * ( 0.2340 * TMath::Gaus(E,4620.5,sigmaTh)
                              +0.763 * TMath::Gaus(E,4687.0,sigmaTh) )
                     +par[2] * ( 0.769 * TMath::Gaus(E,5804.77,sigmaCm)
                              +0.231 * TMath::Gaus(E,5762.16,sigmaCm ) )
                     +par[4];
  return return_f;
}

// ============== PeakHunt() ================ //
// Peak search in TH1 *hist by using TSpectrum;
// Return a vector of x-axis position of 3 peaks with highest y-values;
// TH1 *hist needs to zoom in a reasonable range to skip noise at low ADC channel. 
std::vector<Double_t> PeakHunt(TH1 *hist){
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
  // extract the first 3 elements (3 highest peak) and put their x-values in a new vector top_xpeaks;
  // then order top_xpeaks from lowest to highest based on their values;
  // extract the first and last bin, which should be centers of the first and last true alpha peaks;
  std::sort(peaks.begin(), peaks.end(), [](const auto &m, const auto &n){return m.second > n.second;});
  std::vector<Double_t> top_xpeaks;
  if(npeaks>=3){
    for(int ipeak=0;ipeak<3;ipeak++){
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

// ============== tasf() =================== //
// Return a well defined TF1 for TH1 *h fit
// TH1 *h needs to zoom in a reasonable range to skip noise at low ADC channel. 
// max and min should be centroids of the first and the last true alpha peaks from the source.
// if max or min <0, which means PeakHunt() never been called. Call PeakHunt() to obtain max and min values.
// Options:
// c: open 2 linear calibration parameters. These two parameters are fixed as default.
// l: set the TF1 math formula to LowE_Fun. HighE_Fun as default.
TF1 *tasf(TH1 *h, const char* name = "tas", Double_t min=-1, Double_t max=-1, Option_t *opt=""){
 
  int xbinfirst = h->GetXaxis()->GetFirst();
  int xbinlast  = h->GetXaxis()->GetLast();
  double fit_xlower = h->GetBinLowEdge(xbinfirst); 
  double fit_xupper = h->GetBinLowEdge(xbinlast) + h->GetBinWidth(xbinlast); 
  double ymax = h->GetMaximum(); 

  TString sopt(opt);
  sopt.ToLower();
  sopt.ReplaceAll(" ","");
  
  TF1 *fx = 0;
  Int_t Npx = Int_t((xbinlast-xbinfirst+1)*10);
  
  // Choose the fitting formula
  if(sopt.Contains("l")){
    fx = new TF1(name, TripleAlphaLowE_Fun, fit_xlower, fit_xupper, 9);
    fx->SetParNames("Gd","Th","Cm","fwhmCm","bg","Gd_n","Th_n","offset","gain");
    fx->SetParLimits(5,0.8,1.3);
    fx->SetParLimits(6,0.8,1.2);
    fx->SetNpx(Npx);
    fx->SetParameters(50,50,50,25,0,1,1,0,20);
    fx->SetParLimits(0,0,ymax*10);
    fx->SetParLimits(1,0,ymax*10);
    fx->SetParLimits(2,0,ymax*10);
    fx->SetParLimits(3,20,200);
    fx->SetParLimits(4,0,ymax);    // bg should not be higher than a true peak.
  }else{ // default fit with TripleAlphaHighE_Fun
    fx = new TF1(name, TripleAlphaHighE_Fun, fit_xlower,fit_xupper, 9);
    fx->SetParNames("Pu","Am","Cm","fwhmCm","bg","Pu_n","Am_n","offset","gain");
    fx->SetParLimits(5,0.8,1.3);
    fx->SetParLimits(6,0.8,1.2);
    fx->SetNpx(Npx);
    fx->SetParameters(50,50,50,25,0,1,1,0,20);
    fx->SetParLimits(0,0,ymax*10);
    fx->SetParLimits(1,0,ymax*10);
    fx->SetParLimits(2,0,ymax*10);
    fx->SetParLimits(3,20,200);
    fx->SetParLimits(4,0,ymax);    // bg should not be higher than a true peak.
  }
  
  // Need calibration for the current histogram or not
  if(sopt.Contains("c")){
    // TODO: chekc the size of peaks, if peaks.size() == 0, which means the zoom in range is wrong
    // less than 3 peaks hunting in the hist. Reset the range! 
    if(max<0 || min<0){
      std::vector<double> xpeaks = PeakHunt(h);
      min = xpeaks.front();
      max = xpeaks.back();
      //printf("Inside loop:::max = %f\t min = %f\n", max, min);
    }
    //printf("Outside loop:::max = %f\t min = %f\n", max, min);
    Double_t gain, offset;
    if(sopt.Contains("l")){
      gain = (5804.77-3182.69)/(max-min);
    }else{
      gain = (5804.77-5156.59)/(max-min);
    }
    offset = 5804.77 - gain*max;
    fx->SetParameter(7,offset);  
    fx->SetParameter(8,gain); 
  }else{
    fx->FixParameter(7,0);  // Calibration offset = 0
    fx->FixParameter(8,1);  // Calibration gain   = 1
  }
  return fx;
}


// ================================ After this, need GRSISort Structure ======================== //
// Initialize global varaibles
void Initialize(){
  hlist = new TList; 
  flist = new TList; 
  sumc = new TH2D("sumc", "Channel vs Uncalibarted Charge", 1e4,0,1e4, 1e3,0,1e3); 
  sume = new TH2D("sume", "Channel vs Calibarted Energy"  , 1e4,0,1e4, 1e3,0,1e3); 
}


// Make Hist from Fragment.root file for channels input in main() only
// Hist is for uncalibrated charge
// Save hists into the global TList *hlist;
void MakeHist(std::string infile, const char* calfile, int minCH, int maxCH){

  TChain *fragtree = new TChain("FragmentTree");
  int size = infile.find_last_of("/");
  std::string directory = infile.substr(0,size);
  int num = infile.find_last_of("_");
  std::string runnumber = infile.substr(num-5,5);
  std::cout << "Directory: " << directory << ", Run number: " << runnumber << std::endl;

  //Finds all subruns for a specific run
  std::vector<std::string> runlist;
  DIR * pDIR;
  struct dirent * entry;
  if ((pDIR = opendir(directory.c_str()))) {
    while ((entry = readdir(pDIR))) {
      if (strstr(entry->d_name, runnumber.c_str())) {
        std::string file = directory;
  file.append("/");
  file.append(entry->d_name);
  if(strstr(file.c_str(),"fragment")) runlist.push_back(file);
      }
    }
    closedir(pDIR);
  }   
  std::sort(runlist.begin(),runlist.end()); // Puts subruns in order
  for(int i = 0; i<runlist.size(); i++) {
    fragtree->Add(runlist.at(i).c_str());
  }

  long nentries = fragtree->GetEntries();
  TFragment * tfrag = 0;
  fragtree->SetBranchAddress("TFragment", & tfrag);
  TChannel::ReadCalFile(calfile);

  long xentry = 0;
  for(xentry;xentry<nentries;xentry++){
    fragtree->GetEntry(xentry);
    int chnum = tfrag->GetChannelNumber();
    if(chnum>=minCH && chnum<=maxCH){
      double charge = tfrag->GetCharge();
      if(calmap.empty()){
        sumc->Fill(charge,chnum);
      }else{
        double energy = charge*calmap[chnum].first + calmap[chnum].second;
        sume->Fill(energy,chnum);
      }
    }
    if((xentry%10000)==0){
      printf("Making Hist on entry: %lu / %lu \r", xentry, nentries);
      fflush(stdout);
    }
  }
  printf("Making Hist DONE!  Entry: %lu / %lu \n", xentry, nentries);

  if(calmap.empty()){
    for(int i=0;i<1000;i++){
      if(i>=minCH && i<=maxCH) {
        TH1D *hist = sumc->ProjectionX(Form("Charge_CH%i",i),i+1,i+1);
        hlist->Add(hist);
      }
    }
  }
}

// Run this after "MakeHist()"
// Fit hists generated in MakeHist()
// Fill hists and TF1*fit to Hist.root
// Save gain[] and offsets[] array to Calibration.txt
void CalHist(int minCH, int maxCH, Option_t *opt=""){

  std::vector<int> vec_chan; 
  std::vector<double> vec_gain; 
  std::vector<double> vec_offs; 
  std::vector<double> vec_rePu; 
  std::vector<double> vec_reAm; 
  std::vector<double> vec_reCm; 

  for(int ich=minCH; ich<=maxCH; ich++){
    TH1D *hist = (TH1D *)hlist->FindObject(Form("Charge_CH%i",ich));
    vec_chan.push_back(ich);
    std::vector<Double_t> top_xpeaks = PeakHunt(hist);
    //for(int i=0;i<top_xpeaks.size();i++){
    //  printf("HELLO WORLD\txpeak[%i] = %f\n",i, top_xpeaks.at(i));
    //}
    if(top_xpeaks.size()<3){
      //found # of peaks < 3. something wrong with the current hist
      vec_gain.push_back(1);
      vec_offs.push_back(0);
      vec_reCm.push_back(-1);
      vec_reAm.push_back(-1);
      vec_rePu.push_back(-1);
    }else{
      Double_t min = top_xpeaks.front();
      Double_t max = top_xpeaks.back();
      double xwidth = (max-min)/2.;
      hist->GetXaxis()->SetRangeUser(min-xwidth, max+xwidth);  
      TF1 *fc = tasf(hist, Form("fc_CH%i",ich), min,max,"c");
      hist->Fit(fc,"LQ");     
      hist->Fit(fc,"LQ");     
      flist->Add(fc);
      double gain = fc->GetParameter("gain");
      double offset = fc->GetParameter("offset");
      double reCm = fc->GetParameter("fwhmCm");
      double reAm = reCm * fc->GetParameter(6);
      double rePu = reCm * fc->GetParameter(5);
      vec_gain.push_back(gain); 
      vec_offs.push_back(offset); 
      vec_reCm.push_back(reCm); 
      vec_reAm.push_back(reAm);
      vec_rePu.push_back(rePu);
    }
    calmap[vec_chan.back()] = std::make_pair(vec_gain.back(), vec_offs.back());
  } // hist loop end 

  // Output With
  std::cout << "CHANNEL" << "\t" 
            << "FWHM(Pu)"<< "\t\t"
            << "FWHM(Am)"<< "\t\t"
            << "FWHM(Cm)"<< "\t\t"
            << "GAIN"    << "\t\t"
            << "OFFSET"  << std::endl;
  // Write parameters into the txt file
  std::ofstream outfile;
  outfile.open("Calibration.txt");
  outfile << "float GAIN[" << vec_gain.size() <<"] = {" << vec_gain[0];
  for(int i=0;i<vec_chan.size();i++){
    std::cout << vec_chan[i] << "\t"
              << vec_rePu[i] << "\t\t"
              << vec_reAm[i] << "\t\t"
              << vec_reCm[i] << "\t\t"
              << vec_gain[i] << "\t\t"
              << vec_offs[i] << std::endl;
    if(i>0){
      outfile << ", " << vec_gain[i];
    }
  }
  outfile << "};\n" << "float OFFSET[" << vec_offs.size() << "] = {" << vec_offs[0];
  for(int i=1;i<vec_chan.size();i++){
    outfile << ", " << vec_offs[i];
  }
  outfile << "};";
  outfile.close();
  std::cout << "Writing gains and offsets to " << "Calibration.txt" << std::endl;

}


// =============== mian() ================ //
// Input File:
// 1. FragmentTree.root.file
// 2. Calibration.file
// 3. Min channel number
// 4. Max channel number
int main(int agrc, char **agrv){
  char const *infile;
  char const *calfile;
  char const *outfile;
  char const *ChMin;
  char const *ChMax;
  int ChMin_int;
  int ChMax_int;


  // Finds the GRSISYS path to be used by other parts of the grsisort code
  std::string grsi_path = getenv("GRSISYS"); 
  if(grsi_path.length() > 0) {
    grsi_path += "/"; 
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);

  if(agrc<5){
    std::cout << "Insufficient arguments, provide fragmenttree file + CalibrationFile + min_Ch + max Ch" << std::endl;
    return 0;
  }else if(agrc==5){
    infile  = agrv[1];
    calfile = agrv[2];
    ChMin   = agrv[3];
    ChMax   = agrv[4];
  }else{
    printf("Too many arguments\n");
    return 0;
  }

  TParserLibrary::Get()->Load();
  
  ChMin_int = atoi(ChMin);
  ChMax_int = atoi(ChMax);
 
  Initialize(); 
  MakeHist(infile,calfile,ChMin_int,ChMax_int);
  CalHist(ChMin_int,ChMax_int);
  MakeHist(infile,calfile,ChMin_int,ChMax_int);
  
  outfile = "Hist.root";
  TFile *newf = new TFile(outfile, "recreate");
  newf->cd();
  sumc->Write();
  sume->Write();
  hlist->Write();
  flist->Write();
  newf->Close();
  printf("Input file:%s\nCalibration file: %s\nOutput File: %s\nStarting Channel: %s\nEnding Channel: %s\n",infile,calfile,outfile,ChMin,ChMax);
  

  return 0;
}














