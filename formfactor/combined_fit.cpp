#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "formfactor.cpp"
#include "k+k-.cpp"
#include "k0k0.cpp"

// definition of shared parameter charged mode
const int n0 = 10;
int ipar[n0] = { 0, 1, 2, 3, 4, 5 ,6, 7, 8, 9 };

struct GlobalChi2 {
   const  ROOT::Math::IMultiGenFunction *fChi2_1;
   const  ROOT::Math::IMultiGenFunction *fChi2_2;
   const  ROOT::Math::IMultiGenFunction *fChi2_3;

   GlobalChi2(  ROOT::Math::IMultiGenFunction & f1, ROOT::Math::IMultiGenFunction & f2, ROOT::Math::IMultiGenFunction & f3) : 
                                         fChi2_1(&f1), fChi2_2(&f2), fChi2_3(&f3) {}

   double operator() (const double *par) const {
      double p[n0];
      for (int i = 0; i < n0; ++i) 
          p[i] = par[ ipar[i] ];
      return (*fChi2_1)(p) + (*fChi2_2)(p) + (*fChi2_3)(p);
   }
};

ROOT::Fit::FitResult combinedFit(double end) { //end - граница (в МэВ), до которой будет производиться фит

  TGraphAsymmErrors* h0 = fileKK();
  TGraphAsymmErrors* h1 = getGraph(0);
  TGraphAsymmErrors* h2 = getK0K0Phi();

  TF1* f0 = MDVM::Cross_Section(1);
  TF1* f1 = MDVM::Cross_Section(0);
  TF1* f2 = MDVM::Cross_Section(0);

  // perform now global fit
  ROOT::Math::WrappedMultiTF1 wf0(*f0,1);
  ROOT::Math::WrappedMultiTF1 wf1(*f1,1);
  ROOT::Math::WrappedMultiTF1 wf2(*f2,1);

  ROOT::Fit::DataOptions opt;
  
  ROOT::Fit::DataRange range0;
  ROOT::Fit::DataRange range1;
  ROOT::Fit::DataRange range2;
  
  // set the data range
  range0.SetRange(1.043, end);
  range1.SetRange(1.1, end);
  range2.SetRange(1.0, 1.1);
  
  ROOT::Fit::BinData data0(opt,range0);
  ROOT::Fit::BinData data1(opt,range1);
  ROOT::Fit::BinData data2(opt,range2);
  
  ROOT::Fit::FillData(data0, h0);
  ROOT::Fit::FillData(data1, h1);
  ROOT::Fit::FillData(data2, h2);


  ROOT::Fit::Chi2Function ch0(data0, wf0);
  ROOT::Fit::Chi2Function ch1(data1, wf1);
  ROOT::Fit::Chi2Function ch2(data2, wf2);

  GlobalChi2 chi2(ch0, ch1, ch2);

  ROOT::Fit::Fitter fitter;

  //double par0[Npar] = { 1.139, 1.467, .999, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
  double par0[n0] = { 1.067, 1.28, 1.038, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
  
  
  // create before the parameter settings in order to fix or set range on them
  fitter.Config().SetParamsSettings(n0,par0);
  
  /*
  fitter.Config().ParSettings(0).Fix();
  fitter.Config().ParSettings(1).SetLimits(-0.5, -1.E-5);
  fitter.Config().ParSettings(0).SetStepSize(5);*/

  fitter.Config().MinimizerOptions().SetPrintLevel(0);
  fitter.Config().SetMinimizer("Minuit2","Migrad");

  // fit FCN function directly
  int size = data0.Size() + data1.Size() + data2.Size();
  fitter.FitFCN(n0, chi2, 0, size, true);
  ROOT::Fit::FitResult result = fitter.Result();
  result.Print(std::cout);

  TCanvas * c1 = new TCanvas("Simfit","Simultaneous fit of two histograms",
                             10,10,900,900);

  c1->Divide(1,3);
  for(int i=1; i<=3; i++){
    c1->cd(i)->SetGrid();
    c1->cd(i)->SetLogy();
  }
  gStyle->SetOptFit(1111);

  c1->cd(1);
  f0->SetFitResult( result, ipar);
  f0->SetRange(range0().first, range0().second);
  f0->SetLineColor(kBlue);
  h0->SetTitle("Charged kaons");
  h0->GetListOfFunctions()->Add(f0);
  h0->Draw("ap");

  c1->cd(2);
  f1->SetFitResult( result, ipar);
  f1->SetRange(range1().first, range1().second);
  f1->SetLineColor(kRed);
  h1->SetTitle("Neutral kaons");
  h1->GetListOfFunctions()->Add(f1);
  h1->Draw("ap");
  
  c1->cd(3);
  f2->SetFitResult( result, ipar);
  f2->SetRange(range2().first, range2().second);
  f2->SetLineColor(kGreen);
  h2->SetTitle("Neutral kaons. Peak");
  h2->GetListOfFunctions()->Add(f2);
  h2->Draw("ap");
  
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
