#include "fit.h"

using namespace std;


Double_t Gauss(Double_t* x, Double_t* par){ //N, mean, sigma, const, binWidth
    Double_t f = par[0]*par[4]*Gaus(x[0], par[1], par[2], kTRUE) + par[3] ;
    return f;
    }

double HandleTree::variance(int k, int n){ //variance is sigma^2
    return ((k+1.)*(k+2.))/((n+2.)*(n+3.)) - pow( (k+1.)/(n+2.), 2);
}

    
int* HandleTree::getTriggers(){ //YES
    int* triggers = new int [3]; // [0] - TF, [1] - CF, [2] - TF&CF
    
    string add = (conditions!="") ? "&" : "";
    
    triggers[0] = chain->GetEntries((conditions + add + "t==0").c_str());
    triggers[1] = chain->GetEntries((conditions + add + "t==1").c_str());
    triggers[2] = chain->GetEntries((conditions + add + "t==2").c_str());
    
    return triggers;    
}

double** HandleTree::triggerEfficiency(){//YES
    double** eff = new double* [3];  // [0] - TE(эффективность триггера), [1] - dTE+, [2] - dTE-
    for( int i=0; i<3; i++)
        eff[i] = new double [3];     // [0] - tf, [1] - cf, [2] - general
    
    int* triggers = this->getTriggers();
    
    int T = triggers[0];
    int C = triggers[1];
    int TC = triggers[2];
    
    eff[0][0] = (TC+1)/double(C+TC+2);
    eff[0][1] = (TC+1)/double(T+TC+2);
    eff[0][2] = 1 - (1 - eff[0][0])*(1 - eff[0][1]);
    
    eff[1][0] = sqrt( variance(TC, C+TC) );
    eff[2][0] = eff[1][0];
    
    eff[1][1] = sqrt( variance(TC, T+TC) );
    eff[2][1] = eff[1][1];
    
    eff[1][2] = pow( (1-eff[0][0])*eff[1][1], 2) + pow( (1-eff[0][1])*eff[1][0], 2);
    eff[2][2] = eff[2][0];
    
    return eff;    
}

double* HandleTree::getRegistrationEfficiency(){ //YES
    double* regEff = new double [3]; // [0] - efficiency, [1] - -dEff, [2] - +dEff;
    
    if( fitF==NULL )
        makeFit();
    
    int N0 = 20000; //!!перенести в аргументы метода
    regEff[0] = double(fitF->GetParameter(0))/N0;
    regEff[1] = double(fitF->GetParError(0))/N0;
    regEff[2] = regEff[1];
    
    
    return regEff;
}

void HandleTree::Merge(string inPath){ //YES
    chain->Add(inPath.c_str());
    delete fitF; // удалить старые значения полей, т.к. информация становится устаревшей
    delete hist;
    fitF = NULL;
    hist = NULL;
    
    return;
}

void HandleTree::makeHist(){ //YES
    if(hist!=NULL){
        delete hist;
        hist = NULL;
    }
   
    hist = new TH1D(Form("hist%.1f",energy), Form("Inv mass %.1f", energy), 50, 450, 550);
    //chain->Draw("m>>hist", conditions.c_str()); //простой метод получения гистограммы из дерева, но плохой (не хочу рисовать)
    
    chain->GetEntries(); //перед GetTree нужно вставить любую операцию работы с чейном, иначе всё рушится (баг рута, скорее всего)
    TTree* t = (TTree*)chain->GetTree()->CopyTree(conditions.c_str());
    // все события в таком дереве подходящие, нужно только пробежаться по ним и записать в гистограмму
   
    TTreeReader theReader(t);
    TTreeReaderValue<Double_t> invariantMass(theReader, "m");

    while(theReader.Next()){ 
        hist->Fill(*invariantMass);   //заполняю гистограмму переменной
    }
    //cout << hist->GetEntries() << endl;    
    return;
}

void HandleTree::makeFit(){ //YES
    
    if(hist==NULL)
        makeHist();
    
    gStyle->TStyle::SetOptFit(1111);
    gStyle->TStyle::SetOptStat(0);
    TF1 line("lineFunc", "pol0", 450., 550.);
    hist->Fit(&line, "MQ");
    
    double binWidth = hist->GetBinWidth(0);
    double xmin = 450;
    double xmax = 550;

    fitF = new TF1(Form("fitF%.1f",energy), Gauss, xmin, xmax, 5);
    fitF->SetParLimits(0, 0.0, 1.E7);
    fitF->SetParLimits(1, xmin, xmax);
    fitF->SetParLimits(2, 0.1, 40.);
    fitF->SetParLimits(3, 0.0, 5.);
    fitF->SetParameters(hist->GetEffectiveEntries(), 497.6, 5., line.GetParameter(0), binWidth);
    fitF->FixParameter(4, binWidth);
    fitF->SetParNames("N", "Mean", "#sigma", "const", "binWidth");
    hist->Fit(fitF, "QLL");
    
    return;
}

