#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <iostream>

#include "formfactor.cpp"
#include "k+k-.cpp"
#include "k0k0.cpp"


using namespace std;

// definition of shared parameter charged mode
const int N = 4;
const int n0 = 10;
int ipar[n0] = { 0, 1, 2, 3, 4, 5 ,6, 7, 8, 9 };

struct GlobalChi2 {
   vector<ROOT::Fit::Chi2Function*> f;

   GlobalChi2(vector<ROOT::Fit::Chi2Function*> ch){
    f.reserve(N);
        for(int i=0; i<N; i++)
            f[i] = ch[i];
   }

   double operator() (const double *par) const {
      double p[n0];
      for (int i = 0; i < n0; ++i) 
          p[i] = par[ ipar[i] ];
      return (*f[0])(p) + (*f[1])(p) + (*f[2])(p);// + (*f[3])(p);
   }
};

ROOT::Fit::FitResult combinedFit(double end) { //end - граница (в МэВ), до которой будет производиться фит

  TGraphAsymmErrors* h[N] = {fileKK(), getGraph(0), getK0K0Phi(), fileKKPeak()};
  TF1* f[N] = {MDVM::Cross_Section(1), MDVM::Cross_Section(0), MDVM::Cross_Section(0), MDVM::Cross_Section(1)};
  
  vector<ROOT::Math::WrappedMultiTF1*> wf(N);
  vector<ROOT::Fit::DataRange> range(N);
  vector<ROOT::Fit::BinData*> data(N);
  vector<ROOT::Fit::Chi2Function*> chi(N);
  ROOT::Fit::DataOptions opt;
  
  // set the data range
  range[0].SetRange(1.01, end+0.03);
  range[1].SetRange(1.05, end);
  range[2].SetRange(1.0, 1.1);
  range[3].SetRange(1.0, 1.1);
  
  for(int i=0; i<N; i++){
    wf[i] = new ROOT::Math::WrappedMultiTF1(*f[i], 1);
    data[i] = new ROOT::Fit::BinData(opt, range[i]);
    chi[i] = new ROOT::Fit::Chi2Function(*data[i], *wf[i]);
    
    ROOT::Fit::FillData(*data[i], h[i]);
  }
  
  GlobalChi2 chi2(chi);
  ROOT::Fit::Fitter fitter;

  double par0[n0] = { 1.067, 1.28, 1.038, -0.025, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
  
  // create before the parameter settings in order to fix or set range on them
  fitter.Config().SetParamsSettings(n0, par0);
  fitter.Config().MinimizerOptions().SetPrintLevel(0);
  fitter.Config().SetMinimizer("Minuit2","Migrad");

  // fit FCN function directly
  int size = 0;
  for(int i=0; i<N; i++)
    size += data[i]->Size();
    
  fitter.FitFCN(n0, chi2, &par0[0], size);
  ROOT::Fit::FitResult result = fitter.Result();
  result.Print(std::cout);

  TCanvas * c1 = new TCanvas("Simfit","Simultaneous fit of 4 histograms", 10, 10, 1400, 900);

  c1->Divide(2,2);
  
  for(int i=1; i<=4; i++){
    c1->cd(i)->SetGrid();
    c1->cd(i)->SetLogy();
  }
  gStyle->SetOptFit(1111);
  
  vector<string> titles = {"Charged Kaons", "Neutral Kaons", "Neutral Kaons. Peak", "Charged Kaons. Peak"};
  
  for(int i=0; i<N; i++){
      c1->cd(i+1);
      f[i]->SetFitResult( result, ipar);
      f[i]->SetRange(range[i]().first, range[i]().second);
      f[i]->SetLineColor(kBlue);
      h[i]->SetTitle(titles[i].c_str());
      h[i]->GetListOfFunctions()->Add(f[i]);
      h[i]->Draw("ap");
      cout << (*chi[i])(result.GetParams()) << endl;
  }
  
  return result;
}


void checker(){
    const int N = 50;
    double e[N];
    double rho[N];
    double omg[N];
    double phi[N];
    ROOT::Fit::FitResult res;
    for(int i=0; i<50; i++){
        e[i] = 1.2 + (2. - 1.2)*i/N;
        res = combinedFit(e[i]);
        rho[i] = res.GetParams()[0];
        omg[i] = res.GetParams()[1];
        phi[i] = res.GetParams()[2];
    }
    TGraph* gr = new TGraph(N, &e[0], &rho[0]);
    TGraph* go = new TGraph(N, &e[0], &omg[0]);
    TGraph* gp = new TGraph(N, &e[0], &phi[0]);
    
    gr->SetLineColor(kBlue);
    go->SetLineColor(kRed);
    gp->SetLineColor(kGreen);
    
    TCanvas* c = new TCanvas("can", "Checker", 800, 500);
    
    gr->Draw();
    go->Draw("same");
    gp->Draw("same");
    return;
}

void writer(ROOT::Fit::FitResult res){
    double e = 0.99;
    int n = 1000;
    ofstream o("fitResult.dat");
    double* par = new double [res.NPar()];
    const double* constpar = res.GetParams();
    for(int i=0; i<res.NPar(); i++)
        par[i] = constpar[i];
    for(int i=0; i<n; i++){
        o << i << '\t' << e*1E3 << '\t' << MDVM::Cross_Section(&e, par, 0) << endl;
        e += (1.3 - 0.99)/n;
    }
        
    for(int i=0; i<=n; i++){
        o << i << '\t' << e*1E3 << '\t' << MDVM::Cross_Section(&e, par, 0) << endl;
        e += (2.1 - 1.3)/n;
    }
    o.close();
    return;
}
