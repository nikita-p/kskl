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
int iparC[10] = { 0, 1, 2, 3, 4, 5 ,6, 7, 8, 9 };

// neutral mode
int iparN[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

struct GlobalChi2 {
   GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,
                ROOT::Math::IMultiGenFunction & f2) :
      fChi2_1(&f1), fChi2_2(&f2) {}

   // parameter vector is first background (in common 1 and 2)
   // and then is signal (only in 2)
   double operator() (const double *par) const {
      double p1[9];
      for (int i = 0; i < 9; ++i) p1[i] = par[iparC[i] ];

      double p2[9];
      for (int i = 0; i < 9; ++i) p2[i] = par[iparN[i] ];

      return (*fChi2_1)(p1) + (*fChi2_2)(p2);
   }

   const  ROOT::Math::IMultiGenFunction * fChi2_1;
   const  ROOT::Math::IMultiGenFunction * fChi2_2;
};

ROOT::Fit::FitResult combinedFit(double end) { //end - граница (в МэВ), до которой будет производиться фит

  TGraphAsymmErrors* hN = getGraph(0);
  TGraphAsymmErrors* hC = fileKK();

  TF1* fC = MDVM::Cross_Section(1);
  TF1* fN = MDVM::Cross_Section(0);

  // perform now global fit
  ROOT::Math::WrappedMultiTF1 wfC(*fC,1);
  ROOT::Math::WrappedMultiTF1 wfN(*fN,1);

  ROOT::Fit::DataOptions opt;
  ROOT::Fit::DataRange rangeC;
  // set the data range
  rangeC.SetRange(1.043, end);
  ROOT::Fit::BinData dataC(opt,rangeC);
  ROOT::Fit::FillData(dataC, hC);

  ROOT::Fit::DataRange rangeN;
  rangeN.SetRange(1., end);
  ROOT::Fit::BinData dataN(opt,rangeN);
  ROOT::Fit::FillData(dataN, hN);

  ROOT::Fit::Chi2Function chi2_C(dataC, wfC);
  ROOT::Fit::Chi2Function chi2_N(dataN, wfN);

  GlobalChi2 globalChi2(chi2_C, chi2_N);

  ROOT::Fit::Fitter fitter;

  const int Npar = 10;
  double par0[Npar] = { 1.139, 1.467, .999, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
  
  // create before the parameter settings in order to fix or set range on them
  fitter.Config().SetParamsSettings(Npar,par0);
  
  /*
  fitter.Config().ParSettings(0).Fix();
  fitter.Config().ParSettings(1).SetLimits(-0.5, -1.E-5);
  fitter.Config().ParSettings(0).SetStepSize(5);*/

  fitter.Config().MinimizerOptions().SetPrintLevel(0);
  fitter.Config().SetMinimizer("Minuit2","Migrad");

  // fit FCN function directly
  fitter.FitFCN(Npar,globalChi2,0,dataC.Size()+dataN.Size(),true);
  ROOT::Fit::FitResult result = fitter.Result();
  result.Print(std::cout);

  TCanvas * c1 = new TCanvas("Simfit","Simultaneous fit of two histograms",
                             10,10,900,900);

  c1->Divide(1,2);
  c1->cd(1)->SetGrid();
  c1->cd(2)->SetGrid();
  c1->cd(1)->SetLogy();
  c1->cd(2)->SetLogy();
  c1->cd(1);  
  gStyle->SetOptFit(1111);

  fC->SetFitResult( result, iparC);
  fC->SetRange(rangeC().first, rangeC().second);
  fC->SetLineColor(kBlue);
  hC->SetTitle("Charged kaons");
  hC->GetListOfFunctions()->Add(fC);
  hC->Draw("ap");

  c1->cd(2);
  fN->SetFitResult( result, iparN);
  fN->SetRange(rangeN().first, rangeN().second);
  fN->SetLineColor(kRed);
  hN->SetTitle("Neutral kaons");
  hN->GetListOfFunctions()->Add(fN);
  
  
  hN->Draw("ap");
  
  
  //Comparison with e+e- --> K0K0 CMD-3 data.
  //Phi K0K0 graph. This data is not fitted.
  //It is needed to expand x and y limits on the bottom graph to see this
  TGraphAsymmErrors* hphik0k0 = phik0k0();
  hphik0k0->Draw("p same"); 
  //
  
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
    int n = 200;
    ofstream o("fitResult.dat");
    double* par;
    const double* constpar = res.GetParams();
    for(int i=0; i<res.NPar(); i++)
        par[i] = constpar[i];
    for(int i=0; i<=n; i++){
        o << i << '\t' << e*1E3 << '\t' << MDVM::Cross_Section(&e, par, 0) << endl;
        e += (2. - 0.99)/n;
    }
    o.close();
    return;
}
