//g++ FitRawHist.cxx -Wl,--no-as-needed `root-config --cflags --libs --glibs` -lSpectrum -lMinuit -lGuiHtml -lTreePlayer -lTMVA -L/opt/local/lib -lX11 -lXpm -O2 -Wl,--copy-dt-needed-entries -L/opt/local/lib -lX11 -lXpm `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -lROOTTPython -o FitRawHist


#include <iostream>  
#include <iomanip>   
#include <fstream>   
#include <cmath>     
#include <algorithm>
#include <string> 
#include <stdio.h>   
#include "TH1.h"     
#include "TF1.h"     
#include "TFile.h"   
#include "TList.h"   
#include "TKey.h"
#include "TSpectrum.h"



TList *hlist;
TList *flist;
std::vector<int> vec_chan; 
std::vector<double> vec_gain; 
std::vector<double> vec_offs; 
std::vector<double> vec_rePu; 
std::vector<double> vec_reAm; 
std::vector<double> vec_reCm;


// =============== Initialize() =================== //
void Initialize(){
  hlist = new TList;
  flist = new TList;
}

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
                     +par[2] * ( 0.769 * TMath::Gaus(E,5804.77,sigmaCm)     // Cm
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
  hist->GetXaxis()->SetRangeUser(50, hist->GetNbinsX()*hist->GetBinWidth(1)-1); // skip the nosiy peak
  Int_t npeaks = s->Search(hist,1,"",0.08); // rough peak-search through the entire hist
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
// d: for double peaks
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
    fx->SetParLimits(3,20,500);
    fx->SetParLimits(4,0,ymax);    // bg should not be higher than a true peak.
  }else{ // default fit with TripleAlphaHighE_Fun
    fx = new TF1(name, TripleAlphaHighE_Fun, fit_xlower,fit_xupper, 9);
    fx->SetParNames("Pu","Am","Cm","fwhmCm","bg","Pu_n","Am_n","offset","gain");
    fx->SetParLimits(5,0.5,1.5);
    fx->SetParLimits(6,0.5,1.5);
    fx->SetNpx(Npx);
    fx->SetParameters(50,50,50,25,0,1,1,0,20);
    fx->SetParLimits(0,0,ymax*10);
    fx->SetParLimits(1,0,ymax*10);                                                                                    
    fx->SetParLimits(2,0,ymax*10);
    fx->SetParLimits(3,20,500);
    fx->SetParLimits(4,0,ymax);    // bg should not be higher than a true peak.
    if(sopt.Contains("d")){ // if there are only two peaks in the spectrum. Assume they are from Am and Cm
      fx->FixParameter(0,0);
      fx->FixParameter(5,1);
    }
  } 
    
  // Need calibration for the current histogram or not
  if(sopt.Contains("c")){
    // TODO: chekc the size of peaks, if peaks.size() == 0, which means the zoom in range is wrong
    // less than 3 peaks hunting in the hist. Reset the range! 
    if(max<0 || min<0){
      std::vector<double> xpeaks = PeakHunt(h);
      min = xpeaks.front();
      max = xpeaks.back();
    }
    Double_t gain, offset;
    if(sopt.Contains("l")){
      gain = (5804.77-3182.69)/(max-min);
    }else{
      if(sopt.Contains("d")){
        gain = (5804.77-5485.56)/(max-min);
      }else{
        gain = (5804.77-5156.59)/(max-min);
      }
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


// =============== ReadHist =================== //
void ReadHist(const char* fname){
  TFile *file = TFile::Open(fname);
  if (!file || file->IsZombie()) {
    printf("Error: cannot open file %s\n", fname);
    return;
  }
  TIter nextKey(file->GetListOfKeys());
  TKey *key = nullptr;    
  while ((key = (TKey *)nextKey())) {
    if (strcmp(key->GetClassName(), "TH1D") != 0) continue;
    TH1D *h = (TH1D *)key->ReadObj();
    if (h->GetEntries()<10) continue;
    hlist->Add(h);
  }
}

// =============== CalRawHist() =================== //
void CalRawHist(){
  TIter nextHist(hlist);
  TH1D *hist = nullptr;  
  
  while((hist = (TH1D *)nextHist())){
    TString hname = hist->GetName();
    TString ch = hname(2,hname.Length()-2);
    vec_chan.push_back(std::stoi(ch.Data()));
    std::vector<Double_t> top_xpeaks = PeakHunt(hist);   
    if(top_xpeaks.size()<2){
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
      TF1 *fc;
      if(top_xpeaks.size()==2){
        fc = tasf(hist, Form("fc_CH%i",std::stoi(ch.Data())), min,max, "cd");
      }else{
        fc = tasf(hist, Form("fc_CH%i",std::stoi(ch.Data())), min,max, "c");
      }
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
      flist->Add(fc);
    }
  } // hist loop over
}

// =============== WriteFile() =================== //
void WriteFile(){
  std::cout << "#CHANNEL" << "\t" 
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


// =============== main() =================== //
// Input File:
// 1. raw_hist.root: made by "RawHistMaker.cxx"
int main(int argc, char **argv){
  
  if(argc<2){
    printf("Input raw_hist.root\ni");
    return 1;
  }
  
  const char *fname = argv[1];
  // Step1: Read all hist;
  Initialize();
  ReadHist(fname); 
  // Step2: Fit each hist
  CalRawHist();
  // Step3: Print fitting results and save gain and offset into .dat file
  WriteFile();
  // Step 4: Write raw histograms + fitting fx into fit_hist.root
  TFile *newf = new TFile("fit_hist.root","recreate");
  newf->cd();
  hlist->Write();
  flist->Write();
  newf->Close();

  return 0;
}
