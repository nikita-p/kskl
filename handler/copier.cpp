#include "copier.h"

std::vector<TGraphAsymmErrors*> Copier::showTriggers(){
    
    int n = vec.size();
    
    
    double** tr;
    double** T = new double* [3];
    double** C = new double* [3];
    double** TC = new double* [3];
    for(int i=0; i<3; i++){
        T[i] = new double [n];
        C[i] = new double [n];
        TC[i] = new double [n];
    }
    double* en = new double [n];
    
    
    for(int i=0; i<n; i++){
        tr = vec[i]->triggerEfficiency();
        en[i] = vec[i]->getEnergy();
        T[0][i] = tr[0][0];   T[1][i] = tr[1][0];  T[2][i] = tr[2][0];
        C[0][i] = tr[0][1];   C[1][i] = tr[1][1];  C[2][i] = tr[2][1];
        TC[0][i] = tr[0][2]; TC[1][i] = tr[1][2]; TC[2][i] = tr[2][2];
    }
    
    TGraphAsymmErrors* gT = new TGraphAsymmErrors(n, en, T[0], 0, 0, T[2], T[1]);
    TGraphAsymmErrors* gC = new TGraphAsymmErrors(n, en, C[0], 0, 0, C[2], C[1]);
    TGraphAsymmErrors* gTC = new TGraphAsymmErrors(n, en, TC[0], 0, 0, TC[2], TC[1]);
    
    gT->SetLineColor(1);
    gC->SetLineColor(2);
    gTC->SetLineColor(3);
    
    gT->SetMarkerColor(1);
    gC->SetMarkerColor(2);
    gTC->SetMarkerColor(3);
    
    
    std::vector<TGraphAsymmErrors*> gr;
    gr.insert( gr.end(), gT);
    gr.insert( gr.end(), gC);
    gr.insert( gr.end(), gTC);
    
    
    return gr;
}


vector<TH1D*> Copier::showHists(){

    vector<TH1D*> hs;
    for(int i=0; i<int(vec.size()); i++){
        vec[i]->makeHist();
        hs.insert(hs.end(), vec[i]->getHist());
    }

    return hs;
}

vector<TF1*> Copier::showFits(){
    
    vector<TF1*> vt;
    for(int i=0; i<int(vec.size()); i++){
        vec[i]->makeFit();
        vt.insert(vt.end(), vec[i]->getFit());
    }
    return vt;
}


