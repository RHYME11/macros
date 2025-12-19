//g++ RawHistMaker.cxx -Wl,--no-as-needed `root-config --cflags --libs --glibs` -lSpectrum -lMinuit -lGuiHtml -lTreePlayer -lTMVA -L/opt/local/lib -lX11 -lXpm -O2 -Wl,--copy-dt-needed-entries -L/opt/local/lib -lX11 -lXpm `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -lROOTTPython -o RawHistMaker
 
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
#include "TS3.h"

TList *hlist;
// ================================ After this, need GRSISort Structure ======================== //
void Initialize(){
  hlist = new TList;
}

// ============================ Make the unclibrated energy ========================================//
// Make uncalibrated histogram
// Analysis TTree
void MakeRawHist(TChain *chain, char const *calfile){
  long nentries = chain->GetEntries();
  TS3 *s3 = NULL;
  chain->SetBranchAddress("TS3", &s3);
  if(chain->FindBranch("TS3")){
    chain->SetBranchAddress("TS3", &s3);
  }else{
    std::cout << "Branch 'TS3' not found! TS3 variable is NULL pointer" << std::endl;
    return;
  }
  
  TH1D *hs[1100]; // # of histograms; we only write non-empty histograms in the TList
  for(int i=0;i<1100;i++){
    hs[i] = new TH1D(Form("hs%i",i),Form("uncalibrated energy histogram at CH %i",i), 4000,0,4000); 
  }
  
  if(TChannel::ReadCalFile(calfile) < 1) {
    std::cout << "No channels found in calibration file " << calfile << "!" << std::endl;
    return;     
  }
  std::cout<<std::endl;
  
  long xentry = 0;
  for(xentry;xentry<nentries;xentry++){                                                                               
    chain->GetEntry(xentry);
    for(int i=0;i<s3->GetSectorMultiplicity();i++){
      TS3Hit *sec_hit = s3->GetSectorHit(i);
      int sec_ch   = sec_hit->GetChannelNumber();
      double sec_c = sec_hit->GetCharge();
      hs[sec_ch]->Fill(sec_c);
    }// i (sector) loop over
    for(int i=0;i<s3->GetRingMultiplicity();i++){
      TS3Hit *ring_hit = s3->GetRingHit(i);
      int ring_ch   = ring_hit->GetChannelNumber();
      double ring_c = ring_hit->GetCharge();
      hs[ring_ch]->Fill(ring_c);
    }// i (ring) loop over
    if((xentry%10000)==0){
      printf("Making Hist on entry: %lu / %lu \r", xentry, nentries);
      fflush(stdout);
    } 
  } // entries loop over 
  for(int i=0;i<1100;i++){
    if(hs[i]->GetEntries()>10){ // histogram must not be empty 
      hlist->Add(hs[i]);
    }
  }
  printf("Making Raw Hist DONE!  Entry: %lu / %lu \n", xentry, nentries);
}


// ====================================== main() ==========================================//
// argv1: CalibrationFile
// argv2...: AnalysisTree File Path
int main(int argc, char** argv){
  if(argc<2){
    printf("Input Calibration file and Analysistree file paths");
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

  // Step 2: make uncalibrate energy
  char const *calfile = argv[1];
  Initialize();
  MakeRawHist(chain, calfile);

  // Step 3: Write raw histograms into output.root
  TFile *newf = new TFile("raw_hist.root","recreate");
  newf->cd();
  hlist->Write();
  newf->Close();  

  return 0;
}
