#include <iostream>
#include <TFile.h>

using namespace std;

TGraphAsymmErrors* fileKK(){
    const int n = 56;
    double e[n];
    double s[2][n];
    double stat, syst;
    
    ifstream o( "k+k-.dat" );
    
    int k;
    for(int i=0; i<n; i++){
        o >> e[i] >> s[0][i] >> stat >> syst;
        s[1][i] = sqrt( stat*stat + syst*syst );
        //cout << stat << "\t" << syst << "\t" << s[1][i] << endl;
    }
    o.close();
        
    TGraphAsymmErrors* g = new TGraphAsymmErrors(n, &e[0], &s[0][0], 0, 0, &s[1][0], &s[1][0]);
    /*
    TCanvas* c = new TCanvas("cankk", "canvaskk", 900, 600); 
    g->Draw("ap");
    c->SetGrid();
    c->SetLogy();*/
    //e[2] = ; s[0][2] = ; s[1][2] = ;
    return g;
}

TGraphAsymmErrors* fileKKPeak(){
    const int n = 24;
    double e[2][n];
    double s[2][n];
    
    ifstream o( "k+k-_koz.dat" );
    
    int k;
    for(int i=0; i<n; i++){
        o >> e[0][i] >> e[1][i] >> s[0][i] >> s[1][i];
        e[0][i] *= 1.E-3;
        e[1][i] *= 1.E-3;
        //cout << stat << "\t" << syst << "\t" << s[1][i] << endl;
    }
    o.close();
        
    TGraphAsymmErrors* g = new TGraphAsymmErrors(n, &e[0][0], &s[0][0], &e[1][0], &e[1][0], &s[1][0], &s[1][0]);
    /*
    TCanvas* c = new TCanvas("cankk", "canvaskk", 900, 600); 
    g->Draw("ap");
    c->SetGrid();
    c->SetLogy();*/
    //e[2] = ; s[0][2] = ; s[1][2] = ;
    return g;
}
