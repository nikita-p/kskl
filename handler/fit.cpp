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
    
    int N0 = chain->GetMaximum("soft_ph"); //!!перенести в аргументы метода
    regEff[0] = double(fitF->GetParameter(0))/N0;
    regEff[1] = double(fitF->GetParError(0))/N0;
    regEff[2] = regEff[1];
    
    
    return regEff;
}

void HandleTree::Merge(string inPath, vector<double> lum){ //YES
    chain->Add(inPath.c_str());
    //energies.insert(energies.end(), atof(inPath.substr( inPath.find_last_of('/') + 2, inPath.length() - 5).c_str()) );
    energies.insert(energies.end(), chain->GetMaximum("be") );
    vector<double> zeros = {0,0,0};
    lums.insert(lums.end(), lum[0]>0 ? lum : zeros);
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
   
    hist = new TH1D(Form("hist%.1f",getEnergy()[0]), Form("Inv mass %.1f", getEnergy()[0]), 50, 450, 550);
    TCanvas c("C", "Can", 500, 400);
    c.SetBatch(kTRUE);
    chain->Draw(Form("m>>hist%.1f", getEnergy()[0]), conditions.c_str()); //простой метод получения гистограммы из дерева, но плохой (не хочу рисовать)
    c.SetBatch(kFALSE);
    gROOT->SetBatch(kFALSE);
    /*
    chain->GetEntries(); //перед GetTree нужно вставить любую операцию работы с чейном, иначе всё рушится (баг рута, скорее всего)
    TTree* t = (TTree*)( chain->GetTree()->CopyTree(conditions.c_str()) );
    // все события в таком дереве подходящие, нужно только пробежаться по ним и записать в гистограмму
   
    TTreeReader theReader(t);
    TTreeReaderValue<Double_t> invariantMass(theReader, "m");

    int counter = 0;
    while(theReader.Next()){ 
        hist->Fill(*invariantMass);   //заполняю гистограмму переменной
        counter++;
    }
    cout << counter << endl;    */
    return;
}

void HandleTree::makeFit(){ //YES
    
    makeHist();
    
    gStyle->TStyle::SetOptFit(1111);
    gStyle->TStyle::SetOptStat(0);
    TF1 line("lineFunc", "pol0", 450., 550.);
    hist->Fit(&line, "MQ");
    
    double binWidth = hist->GetBinWidth(0);
    double xmin = 450;
    double xmax = 550;

    fitF = new TF1(Form("fitF%.1f",getEnergy()[0]), Gauss, xmin, xmax, 5);
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

vector<double> HandleTree::getEnergy(){
    double EL = 0;
    double L = 0;
    int n = energies.size();
    for(int i=0; i<n; i++){
        EL += energies[i]*lums[i][0];
        L += lums[i][0];
    }
    vector<double> v(3);
    if( fabs(EL)>0.001 ){
        v[0] = EL/L;
        v[1] = (*max_element(energies.begin(),energies.end()) - EL/L);
        v[2] = (EL/L - *min_element(energies.begin(),energies.end()));
        return v;
    }
    return {-energies[n/2],0,0}; //отрицательная энергия указывает на то, что не все светимости определены корректно
}

vector<double> HandleTree::getLum(){
    int n = lums.size();
    vector<double> L = {0, 0, 0};
    for(int i=0; i<n; i++){
        L[0] += lums[i][0];
        L[1] += lums[i][1];
        L[2] += lums[i][2];
        }
    return L;
}
