#include "fit.h"
#include "copier.h"

void gogogo(){
    
    Copier a("Outputs/trees/2011/");
    vector<TH1D*> tt = a.showHists();

    TCanvas* c = new TCanvas("cvs", "CAAANVAAAS", 1200, 700);
    c->Divide(8,5);
    for(int i=0; i<int(tt.size()); i++){
        c->cd(i+1);
        tt[i]->Draw();
        c->Update();
    }
    
    return;
}
