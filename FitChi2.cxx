//g++ FitChi2.cxx -o bin/FitChi2 `root-config --cflags --libs`
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cmath>

#include <TGraph.h>
#include <TF1.h>
#include <TList.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h>

// ============= new data structure ============= //
struct Chi2Data{
  int DOE; // degree of freedom
  std::vector<double> taus;
  std::vector<double> chi2s;
};





// ============= Read File ====================//
Chi2Data ReadChi2File(const std::string& filename){
  
  Chi2Data result;

  std::ifstream file(filename);
  if(!file.is_open()){
    std::cerr << "Error opening file: " << filename << std::endl;
    return Chi2Data{};
  }
  double lowerE = 0;
  double upperE = 0;
  double binwidth = 0;
  std::vector<double> taus;
  std::vector<double> chi2s;
    
  std::string line;
  bool header_parsed = false;
  
  while (std::getline(file, line)){
    // skip empty space
    if(line.empty()) continue;

    // readout "lowwrE", "upperE" and "binwidth"
    if (!header_parsed && line[0] == '#') {
      std::stringstream ss(line);
      std::string tmp;
      while (ss >> tmp) {
        if (tmp.find("lowerE=") != std::string::npos) {
          lowerE = std::stod(tmp.substr(tmp.find('=') + 1));
        } else if (tmp.find("upperE=") != std::string::npos) {
          upperE = std::stod(tmp.substr(tmp.find('=') + 1));
        } else if (tmp.find("binwidth=") != std::string::npos) {
          std::string bw = tmp.substr(tmp.find('=') + 1);
          size_t pos = bw.find_first_not_of("0123456789.");
          binwidth = std::stod(bw.substr(0, pos));
        }
      }
      header_parsed = true;
      continue;
    }
    
    result.DOE = (int)(upperE - lowerE)/binwidth;

    // skip the title line
    if (line.find("tau") != std::string::npos && line.find("chi2") != std::string::npos)
      continue;
    
    // readout tau-values and chi2-values
    std::stringstream ss(line);
    double tau, chi2, scaling;
    if (ss >> tau >> chi2 >> scaling) {
      result.taus.push_back(tau);
      result.chi2s.push_back(chi2);
    }
  } // while loop over
    
  return result;
}


// =========== Return tau corresponding min_chi2 (stationary point) ============== //
double FindTau_SP(TF1 *fx){
  double c = fx->GetParameter(1);  
  double b = fx->GetParameter(2); 
  double a = fx->GetParameter(3); 

  double x = sqrt((b*b - 3*a*c)/(9*a*a)) - b/(3*a); // d(pol3)/dx = 0 ==> x = formular
  return x;
}


// =========== Return tau corresponding a given chi2 (root) ============== //
std::vector<double> FindTau(TF1 *fx, double chi2){
  
  double d = fx->GetParameter(0); 
  double c = fx->GetParameter(1); 
  double b = fx->GetParameter(2); 
  double a = fx->GetParameter(3); 
  d = d-chi2;

  double alpha = (b*c)/(6*a*a) - (b*b*b)/(27*a*a*a) - d/(2*a);
  double beta  = c/(3*a) - (b*b)/(9*a*a);
  double arcx = alpha/pow(-beta,1.5);
  
  double x1 = -b/(3*a) + 2*sqrt(-beta)*cos(acos(arcx)/3.);
  double x2 = -b/(3*a) + 2*sqrt(-beta)*cos((acos(arcx)+2*3.1415926)/3.);
  double x3 = -b/(3*a) + 2*sqrt(-beta)*cos((acos(arcx)-2*3.1415926)/3.);

  std::vector<double> result;
  if(x1>fx->GetXmin() && x1<fx->GetXmax()) result.push_back(x1);
  if(x2>fx->GetXmin() && x2<fx->GetXmax()) result.push_back(x2);
  if(x3>fx->GetXmin() && x3<fx->GetXmax()) result.push_back(x3);

  std::sort(result.begin(), result.end());
  return result; // result[0] < result[1] always
}


// ========= main function ============= //
int main(int argc, char **argv){
  
  std::vector<std::string> filenames;

  if(argc<2){
    printf("Add Inputs!\n");
    return 1;
  }else{
    for(int i=1;i<argc;i++){
      filenames.push_back(argv[i]);
    } 
  }

  
  std::ofstream outfile("fitting_results.txt");
  TList *glist = new TList;
  TList *flist = new TList;
  
  for(size_t i=0;i<filenames.size();i++){
    // ==== Read Chi2.txt ==== //
    const std::string& fname = filenames[i];
    Chi2Data data = ReadChi2File(fname);
    int ndata = data.taus.size();
    int doe = data.DOE; // degree of freedom
    double tau_min = *std::min_element(data.taus.begin(), data.taus.end());
    double tau_max = *std::max_element(data.taus.begin(), data.taus.end());

    // ==== Fit chi2 vs tau (with pol3) ==== //
    TGraph *gr = new TGraph(ndata, &data.taus[0], &data.chi2s[0]);
    gr->SetName(Form("gr%zu",i));
    TF1 *fx = new TF1(Form("fx%zu",i), "pol3", tau_min, tau_max);
    gr->Fit(fx,"Q"); 
    fx->SetTitle(Form("chi2 = %.4f + (%.4f*x) + (%.4f*x^2) + (%.4f*x^3)",
                      fx->GetParameter(0),
                      fx->GetParameter(1),
                      fx->GetParameter(2),
                      fx->GetParameter(3)));
    glist->Add(gr);
    flist->Add(fx);
    
    // ==== numerical calculation ==== // 
    double chi2_min = fx->GetMinimum(tau_min,tau_max);
    double tau = FindTau_SP(fx); // tau corresponding to chi2_min
    // === if chi2_min values are different obtained through two methods, trust my calculation ==== //
    double dif = fabs(chi2_min - fx->Eval(tau))/chi2_min;
    if(dif>0.01){
      printf("Something wrong with the fit%zu! [min, eval] = [%f, %f], tau = %.2f\n",i, chi2_min, fx->Eval(tau), tau);
      chi2_min = fx->Eval(tau);
    }
    std::vector<double> taus = FindTau(fx, chi2_min+1); // 2 tau-values corresponding to chi2_min+1
    if(taus.size() != 2){
      printf("something wrong with FindTau function in fit%zu\n",i);
      continue;
    }
    std::vector<double> sigma;
    std::vector<double> red_sigma;
    for(int j=0;j<2;j++){
      sigma.push_back(fabs(taus[j] - tau));
      red_sigma.push_back(sigma[j] * pow((tau*tau)/doe, 0.5));
    }   
    
    // ==== record numerical results ==== //
    outfile << "# " << fname << std::endl;
    outfile << "[tau, chi2_min] = " << tau << ", " << chi2_min << std::endl;
    outfile << "[tau1, tau2, chi2_min+1] = " << taus[0] << ", " << taus[1] << ", " << chi2_min+1 << std::endl;
    outfile << "[sig1, sig2] = " << sigma[0] << ", " << sigma[1] << std::endl;
    outfile << "[red_sig1, red_sig2] = " << red_sigma[0] << ", " << red_sigma[1] << std::endl;
    outfile << std::endl; 
    
  } // argv loop over
  
  TFile *newf = new TFile("fitting_results.root","recreate");
  glist->Write();
  flist->Write();
  newf->Close();
  
  return 0;
}















