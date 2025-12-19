//g++ src/HistMakers.cxx -Wl,--no-as-needed `root-config --cflags --libs --glibs` -lSpectrum -lMinuit -lGuiHtml -lTreePlayer -lTMVA -L/opt/local/lib -lX11 -lXpm -O2 -Wl,--copy-dt-needed-entries -L/opt/local/lib -lX11 -lXpm `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -lROOTTPython -o HistMakers


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
#include "TS3.h"
#include "TRandom.h"
#include "TRandom3.h"






TList *hlist;                      
TRandom3 rand3;
// ================================ After this, need GRSISort Structure ======================== //
void Initialize(){                 
  hlist = new TList;               
} 


// ========================= GetS3Position() =======================//
TVector3 GetS3Position(int det, int ring, int sec, bool rings_facing_target = false, bool smear = false) {
  
  if(det < 1 || det > 2 || ring < 0 || ring > 23 || sec < 0 || sec > 31)
    return TVector3(0.0,0.0,0.0);
   
  TVector3 pos(1.0,1.0,0.0); //X and Y arbitrary, Z must be zero
  double rOff = 0;
  if(smear)
  rOff = rand3.Rndm() - 0.5;
  pos.SetPerp(ring + 11.5 + rOff);
    
  double pOff = 0; 
  if(smear)
    pOff = rand3.Rndm() * 11.5 - 11.5/2.;
 
  //JW: sector angles determined based on the source rotation test
  //the 9 sector offset is (I think) related to the 90 degree 
  //difference between the EMMA and Bambino chamber mounts
  //these angles are defined for detector 1 (downstream)
  //with x axis +ve in the direction of TIGRESS position 5
  pos.SetPhi((sec+9.0)*TMath::TwoPi()/32.0);
 
  //S3 detectors have opposite orientation, so flip detector 2
  if(det == 2){
    pos.RotateY(TMath::Pi());
  }
 
  //rotate both detectors based on the Bambino chamber rotation
  pos.RotateZ(-22.5*TMath::DegToRad() + pOff*TMath::DegToRad());
 
  if(rings_facing_target)
    pos.RotateY(TMath::Pi());
 
  pos.SetZ(33.0);
 
  return pos;
  
}
// ============================ Make the unclibrated energy ========================================//
// Make uncalibrated histogram    
// Analysis TTree                 
void MakeHist(TChain *chain, char const *calfile){
  long nentries = chain->GetEntries();
  TS3 *s3 = NULL;                 
  chain->SetBranchAddress("TS3", &s3);
  if(chain->FindBranch("TS3")){   
    chain->SetBranchAddress("TS3", &s3);
  }else{                          
    std::cout << "Branch 'TS3' not found! TS3 variable is NULL pointer" << std::endl;
    return;                       
  }                               
    
  // ~~~~~~~~~~~~~~~~ Hists Definetion ~~~~~~~~~~~~~~~~~~~~~~~ //
  TH2D *s3_XY[2];
  TH1D *dthist[2];
  TH2D *s3_XY_sec[2][32]; 
  for(int i=0;i<2;i++){
    s3_XY[i] = new TH2D(Form("s3_XY_Det%i",i),Form("s3 XY Det%i",i), 250, -40.0, 40.0, 250,-40.0,40.0);
    dthist[i] = new TH1D(Form("dthist_Det%i",i),Form("dt(ns) = Ring.T - Sec.T in Det%i",i), 6000,-3000,3000);
    for(int j=0;j<32;j++){
      s3_XY_sec[i][j] = new TH2D(Form("s3_XY_Det%i_sec%i",i,j),Form("s3 XY at Secctor%i Det%i",j,i), 250, -40.0, 40.0, 250,-40.0,40.0);
    }
  }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
                                  
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
      int sec_det  = sec_hit->GetDetector(); 
      int sec      = sec_hit->GetSector();
      double sec_c = sec_hit->GetCharge();
      double sec_t = sec_hit->GetTime();
      if(sec_c<50) continue;
      for(int j=0;j<s3->GetRingMultiplicity();j++){
        TS3Hit *ring_hit = s3->GetRingHit(j);
        int ring_det  = ring_hit->GetDetector();
        int ring      = ring_hit->GetRing();
        double ring_c = ring_hit->GetCharge();
        double ring_t = ring_hit->GetTime();
        if(ring_c<50) continue;
        if(sec_det != ring_det) continue;
        double dt = ring_t - sec_t;
        dthist[sec_det-1]->Fill(dt);
        TVector3 posS3Smear = GetS3Position(sec_det,ring,sec,false,true);
        posS3Smear.SetX(posS3Smear.X() + 1.35);
        posS3Smear.SetY(posS3Smear.Y() - 0.65);
        s3_XY[sec_det-1]->Fill(posS3Smear.X(), posS3Smear.Y());
        s3_XY_sec[sec_det-1][sec]->Fill(posS3Smear.X(), posS3Smear.Y());
      }// j (ring) loop over
    }// i (sector) loop over       
    if((xentry%10000)==0){         
      printf("Making Hist on entry: %lu / %lu \r", xentry, nentries);
      fflush(stdout);              
    }                              
  } // entries loop over           
  for(int i=0;i<2;i++){
    hlist->Add(s3_XY[i]);
    hlist->Add(dthist[i]);
    for(int j=0;j<32;j++){
      hlist->Add(s3_XY_sec[i][j]);
    }
  }
  
  printf("Making Raw Hist DONE!  Entry: %lu / %lu \n", xentry, nentries);
}

// ====================================== main() ==========================================//
// argv1: CalibrationFile
// argv2...: AnalysisTree File Path
int main(int argc, char** argv){
  if(argc<3){
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
 
  // Step 2: make histograms
  char const *calfile = argv[1];
  Initialize();
  MakeHist(chain, calfile);
 
  // Step 3: Write raw histograms into output.root
  TFile *newf = new TFile("hist.root","recreate");
  newf->cd();
  hlist->Write();
  newf->Close();  
 
  return 0;
}
