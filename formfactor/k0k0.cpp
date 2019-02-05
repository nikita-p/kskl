void reader(TTree* t, int n0, double* cs0, double* cs1, double* cs2, double* e0, double* e1, double* e2, double* rc, bool vis){
    
    TTreeReader theReader(t);
    TTreeReaderArray<double> CS(theReader, "cs");
    TTreeReaderArray<double> E(theReader, "E");
    TTreeReaderValue<double> RC(theReader, "radcor");

    int i = n0;
    while(theReader.Next()){
        cs0[i] = CS[0];
        cs1[i] = CS[1];
        cs2[i] = CS[2];
        if(vis==true){
            cs0[i] *= (*RC);
            cs1[i] *= (*RC);
            cs2[i] *= (*RC);
        }
        e0[i] = E[0]*0.002;
        e1[i] = E[1]*0.002;
        e2[i] = E[2]*0.002;
        rc[i] = *RC;
        i++;
    }

    return;
}

TGraphAsymmErrors* getGraph(bool vis){
    
    TFile* f11 = TFile::Open("res11.root");
    TFile* f12 = TFile::Open("res12.root");
    TFile* f17 = TFile::Open("res17.root");
    
    TTree* t11 = (TTree*)f11->Get("t");
    TTree* t12 = (TTree*)f12->Get("t");
    TTree* t17 = (TTree*)f17->Get("t");
    
    int n11 = t11->GetEntries();
    int n12 = t12->GetEntries();
    int n17 = t17->GetEntries();
    
    const int n = n11 + n12 + n17;
    
    double cs0[n];
    double cs1[n];
    double cs2[n];
    double e0[n];
    double e1[n];
    double e2[n];
    double rc[n];
    
    reader(t11, 0, &cs0[0], &cs1[0], &cs2[0], &e0[0], &e1[0], &e2[0], &rc[0], vis);
    reader(t12, n11, &cs0[0], &cs1[0], &cs2[0], &e0[0], &e1[0], &e2[0], &rc[0], vis);
    reader(t17, n11+n12, &cs0[0], &cs1[0], &cs2[0], &e0[0], &e1[0], &e2[0], &rc[0], vis);
    
    f11->Close();
    f12->Close();
    f17->Close();
    
    TGraphAsymmErrors* gr = new TGraphAsymmErrors(n, &e0[0], &cs0[0], &e2[0], &e1[0], &cs2[0], &cs1[0]);
    /*
    TCanvas* c = new TCanvas("can", "CrossSection", 800, 500);
    gr->Draw("ap");
    c->SetLogy();
    c->SetGrid();
    gStyle->SetOptFit(11111);
    gr->SetTitle("");
    gr->GetXaxis()->SetTitle("#sqrt{s}, GeV");
    gr->GetYaxis()->SetTitle("#sigma, nb");*/
    
    return gr;
}

void writeToFile(TF1* f){
    int n = 2000;
    int n0 = 500;
    int n1 = n - n0;
    
    double emin = 497.6*2;
    double e1 = 1077;
    double emax = 2010;
    double x;
    ofstream o("cs.dat");
    for(int i=0; i<n0; i++){
        x = emin + (e1-emin)*i/n0;
        o << i << '\t' << x << '\t' << f->Eval(x*1e-3) << endl;
    }
    for(int i=0; i<=n1; i++){
        x = e1 + (emax-e1)*i/n1;
        o << i+n0 << '\t' << x << '\t' << f->Eval(x*1e-3) << endl;
    }
    o.close();
    return;
}
