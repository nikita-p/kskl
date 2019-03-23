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
    return g;
}

TGraphAsymmErrors* filePiFormfactor(){
    const int n = 60;
    double e[n];
    double s[3][n];
    
    ifstream o( "pi+pi-_formfactor.dat" );
    
    int k;
    for( int i=0; i<n; i++){
        o >> e[i] >> s[0][i] >> s[1][i] >> s[2][i];
        s[2][i] *= -1; //в этом файле отрицательные нижние ошибки
        e[i] *= 1.E-3;
    }
    o.close();
    
    TGraphAsymmErrors* g = new TGraphAsymmErrors(n, &e[0], &s[0][0], 0, 0, &s[1][0], &s[2][0]);
    return g;
}
