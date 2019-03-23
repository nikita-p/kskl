#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include <iostream>

#include "formfactor.cpp"
#include "k+k-.cpp"
#include "k0k0.cpp"


using namespace std;

// definition of shared parameter charged mode
const int N = 4;
const int n0 = 11;
int ipar[n0] = { 0, 1, 2, 3, 4, 5 ,6, 7, 8, 9, 10 };

void Apply(TGraphAsymmErrors* g, TF1 *f)
 {
    Double_t x,y,exl,exh,eyl,eyh,fxy;
    for(int i=0; i<(g->GetN()); i++) {
       
       g->GetPoint(i,x,y);
       
       exl=g->GetErrorXlow(i);
       exh=g->GetErrorXhigh(i);
       eyl=g->GetErrorYlow(i);
       eyh=g->GetErrorYhigh(i);
       
       fxy = ( y - (f->Eval(x)) );
       g->SetPoint(i,x,fxy);
       g->SetPointError(i,exl,exh,eyl,eyh);
       //g->SetPointError(i,exl,exh,0,0);
       
    }
}

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
      return (*f[0])(p) + (*f[1])(p) + (*f[2])(p) + (*f[3])(p);
   }
};

ROOT::Fit::FitResult combinedFit(double end = 2., bool zeromode = 0) { //end - граница (в МэВ), до которой будет производиться фит
    
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
  range[2].SetRange(1.0, 1.061);
  range[3].SetRange(1.01, 1.061);
  
  for(int i=0; i<N; i++){
    wf[i] = new ROOT::Math::WrappedMultiTF1(*f[i], 1);
    data[i] = new ROOT::Fit::BinData(opt, range[i]);
    chi[i] = new ROOT::Fit::Chi2Function(*data[i], *wf[i]);
    
    ROOT::Fit::FillData(*data[i], h[i]);
  }
  
  GlobalChi2 chi2(chi);
  ROOT::Fit::Fitter fitter;

  double par0[n0] = { 1.067, 1.28, 1.038, -0.0628, -0.0503, -0.19, -0.042, 0.083, -0.077, 0.0, 0.987};
  
  // create before the parameter settings in order to fix or set range on them
  fitter.Config().SetParamsSettings(n0, par0);
  fitter.Config().MinimizerOptions().SetPrintLevel(0);
  fitter.Config().SetMinimizer("Minuit2","Migrad");

  //fitter.Config().ParSettings(0).Fix();
  //fitter.Config().ParSettings(1).Fix();
  //fitter.Config().ParSettings(2).Fix();
  //fitter.Config().ParSettings(3).Fix();
  //fitter.Config().ParSettings(4).Fix();
  //fitter.Config().ParSettings(5).Fix();
  //fitter.Config().ParSettings(6).Fix();
  //fitter.Config().ParSettings(7).Fix();
  //fitter.Config().ParSettings(8).Fix();
  //fitter.Config().ParSettings(9).Fix();
  fitter.Config().ParSettings(10).Fix();
  
  // fit FCN function directly
  int size = 0;
  for(int i=0; i<N; i++)
    size += data[i]->Size();
    
  fitter.FitFCN(n0, chi2, &par0[0], size);
  ROOT::Fit::FitResult result = fitter.Result();
  
  //if(zeromode)
  //  return chi2(result.GetParams());
  
  result.Print(std::cout);

  TCanvas* c[2] = { new TCanvas( "Simfit",  "Simultaneous fit of 4 histograms", 10, 10, 1400, 900), 
                    new TCanvas("Difffit","Differences between fit and points", 10, 10, 1400, 900) };

  c[0]->Divide(2,2);
  c[1]->Divide(2,2);
  TGraphAsymmErrors* hdiff[N];
  for(int i=0; i<N; i++)
    hdiff[i] = (TGraphAsymmErrors*)(h[i]->Clone());
  
  for(int i=1; i<=N; i++){
    c[0]->cd(i)->SetGrid();
    c[1]->cd(i)->SetGrid();
    c[0]->cd(i)->SetLogy();
  }

  gStyle->SetOptFit(1111);
  
  vector<string> titles = {"Charged Kaons", "Neutral Kaons", "Neutral Kaons. Peak", "Charged Kaons. Peak"};
  
  for(int i=0; i<N; i++){
      c[0]->cd(i+1);
      f[i]->SetFitResult( result, ipar);
      f[i]->SetRange(range[i]().first, range[i]().second);
      f[i]->SetLineColor(kBlue);
      h[i]->SetTitle(titles[i].c_str());
      h[i]->GetListOfFunctions()->Add(f[i]);
      h[i]->Draw("ap");
      
      cout << (*chi[i])(result.GetParams()) << endl;
  }
      //cout << chi2(result.GetParams()) << endl;
  //c[0]->Close();
  
  for(int i=0; i<N; i++){
      c[1]->cd(i+1);
      Apply(hdiff[i], f[i]);
      hdiff[i]->SetTitle(titles[i].c_str());
      TLine* l = new TLine(range[i]().first, 0, range[i]().second, 0);
      hdiff[i]->SetMarkerStyle(7);
      hdiff[i]->Draw("ap");
      l->Draw("same");
  }
  c[1]->Close();
  
  void writer(ROOT::Fit::FitResult);
  writer(result);
  return result;//chi2(result.GetParams());
}

/*
void checker(){
    int N[2] = {4, 4};
    pair<double,double> mphi( 1465-25, 1465+25 );
    pair<double,double> wphi(  400-60,  400+60 );
    
    ofstream o("chisquarePhi.dat");
    
    double dmphi = (mphi.second - mphi.first)/(N[0]-1);
    double dwphi = (wphi.second - wphi.first)/(N[1]-1);
    
    double Mphi, Wphi;
    
    for(int i=0; i<N[0]; i++)
      for(int j=0; j<N[1]; j++){
      
        Mphi = mphi.first + i*dmphi;
        Wphi = wphi.first + j*dwphi;
        
        cout << i << '\t' << j << endl;
        o << Mphi << ';' << Wphi << ';' << combinedFit(Mphi, Wphi, 2, true) << endl;
                      
      }
    o.close();
    return;
}
*/
void writer(ROOT::Fit::FitResult res){
    double e = 0.994;
    int n = 1000;
    ofstream o("fitResult.dat");
    double* par = new double [res.NPar()];
    const double* constpar = res.GetParams();
    for(int i=0; i<res.NPar(); i++)
        par[i] = constpar[i];
    for(int i=0; i<n; i++){
        o << i << '\t' << e*1E3 << '\t' << MDVM::Cross_Section(&e, par, 0) << endl;
        e += (1.15 - 0.99)/n;
    }
        
    for(int i=0; i<=n; i++){
        o << i << '\t' << e*1E3 << '\t' << MDVM::Cross_Section(&e, par, 0) << endl;
        e += (2.03 - 1.15)/n;
    }
    o.close();
    return;
}
