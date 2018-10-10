#include "fit.h"
#include "copier.h"

void gogogo(){
    
    Copier a11("Outputs/trees/2011/");
    vector<TH1D*> t = a.showHists();
    vector<TF1*> d = a.showFits();

    TCanvas* c = new TCanvas("cvs", "CAAANVAAAS", 1200, 700);
    c->Divide(8,5);
    for(int i=0; i<int(t.size()); i++){
        c->cd(i+1);
        d[i]->Draw();
        t[i]->Draw("same");
        c->Update();
    }
    
    return;
}
