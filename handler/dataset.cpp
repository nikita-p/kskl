#include "dataset.h"

using namespace std;

bool Dataset::compt(HandleTree* a, HandleTree* b){
    return (a->getEnergy())<(b->getEnergy());
}

bool Dataset::comps(string a, string b){
    double A = atof( a.substr( a.find_last_of('/') + 2, a.length() - 5).c_str() );
    double B = atof( b.substr( b.find_last_of('/') + 2, b.length() - 5).c_str() );
    return A<B;
}

double Dataset::LinearApprox(double E, double x1, double y1, double x2, double y2){
    return y1 + ( (y2 - y1)/(x2 - x1) ) * (E - x1);
}

double* Dataset::LinearApprox(int n, double E, double x1, double* y1, double x2, double* y2){

    double* res = new double [n];
    for(int i=0; i<n; i++)
        res[i] = LinearApprox(E, x1, y1[i], x2, y2[i]);
    return res;
    
}

double Dataset::GetEnergyFromString(string a){
    return atof( a.substr( a.find_last_of('/') + 2, a.length() - 5).c_str() );
}

void Dataset::GetDataVector(){
    void* dir = gSystem->OpenDirectory( dataFolder.c_str() );
    
    if( !dir ){
        std::cout << "Error 1. Bad times" << std::endl;
        return; }
        
    const char* entry;
    string end = ".root";
    TString filename;
        
    entry = gSystem->GetDirEntry(dir);
        
    while( entry ){
        filename = entry;
        if(filename.EndsWith(end.c_str())){
            filename = gSystem->ConcatFileName(dataFolder.c_str(), entry);
            data.insert(data.end(), string(filename));
        }
        entry = gSystem->GetDirEntry(dir);
    }
    sort(data.begin(), data.end(), comps);
    return;
}

void Dataset::GetModelVector(){
    void* dir = gSystem->OpenDirectory( modelFolder.c_str() );
    
    if( !dir ){
        std::cout << "Error 1. Bad times" << std::endl;
        return; }
        
    const char* entry;
    string end = ".root";
    TString filename;
        
    entry = gSystem->GetDirEntry(dir);
        
    while( entry ){
        filename = entry;
        if(filename.EndsWith(end.c_str())){
            filename = gSystem->ConcatFileName(modelFolder.c_str(), entry);
            HandleTree* f = new HandleTree(string(filename), treeName.c_str(), 1); //для моделирования можно любую светимость (главное, одинаковую)
            model.insert(model.end(), f);
        }
        entry = gSystem->GetDirEntry(dir);
    }
    sort(model.begin(), model.end(), compt);
    return;
}

void Dataset::GetLumVector(){
    ifstream flum(lumFile.c_str());
    
    double e;
    while(!flum.eof()){
        double* l = new double [3];
        flum >> e >> l[0] >> l[1];
        //cout << e << '\t' << l[0] << '\t' << l[1] << endl;
        if(flum.eof())
            break;
        l[2] = l[1];
        pair<double, double*> p(e, l);
        lum.insert(lum.end(), p);
    }    
    flum.close();
    return;
}

void Dataset::GetRadcorVector(){
    ifstream frad(radcorFile.c_str());
    
    double e;
    double rad;
    while(!frad.eof()){
        frad >> e >> rad;
        if(frad.eof())
            break;
        pair<double, double> p(e, rad);
        radcor.insert(radcor.end(), p);
    }
    frad.close();
    return;
}

Dataset::Dataset(string model, string data, string lum, string radcor){
    modelFolder = model;
    dataFolder = data;
    lumFile = lum;
    radcorFile = radcor;
    
    GetDataVector();
    GetModelVector();
    GetLumVector();
    GetRadcorVector();
}

double* Dataset::RegistrationEff(double E){
    //gROOT->SetBatch(kTRUE);
    int n = model.size();
    double index[2] = {0, 1};
    double* res = new double [3];
    double x1, x2;
    double* y1;
    double* y2;
    double E0, E1;
    
    for(int i=0; i<n-1; i++){
        E0 = model[i]->getEnergy();
        E1 = model[i+1]->getEnergy();
        if( (E>=E0)&&(E<=E1) ){
            index[0] = i;
            index[1] = i+1;
            break;
        }
    }
    if( E>(model[n-1]->getEnergy()) ){
        index[0] = n-2; index[1] = n-1;
    }
    
    x1 = model[index[0]]->getEnergy();
    x2 = model[index[1]]->getEnergy();
    y1 = model[index[0]]->getRegistrationEfficiency();
    y2 = model[index[1]]->getRegistrationEfficiency();
    res = LinearApprox(3, E, x1, y1, x2, y2);
    //gROOT->SetBatch(kFALSE);
    return res;
}

void Dataset::ClearFits(){
    fits.clear();
    return;
}

void Dataset::ClearHists(){
    hists.clear();
    return;
}

void Dataset::ClearMdata(){
    mdata.clear();
    return;
}

void Dataset::AutoFit(){

    ClearFits();
    ClearHists();
    ClearMdata();
    
    int n = data.size();
    HandleTree* t;
    for(int i=0; i<n; i++){
        t = new HandleTree(data[i], treeName, GetEnergyFromString(data[i]));
        t->makeFit();
        mdata.insert(mdata.end(), t);
        fits.insert(fits.end(), t->getFit());
        hists.insert(hists.end(), t->getHist());
    }
    return;
}

void Dataset::PowerFit(int startPosition = 0){

    ClearFits();
    ClearHists();
    ClearMdata();
    
    int n = data.size();
    HandleTree* t;
    double E;
    bool do_merge = false;
    for(int i=startPosition; i<n; i++){
        E = GetEnergyFromString(data[i]);
        if( fabs(GetLuminosity(E))<0.0001 ) //нет светимости для этого файла -> не юзать его совсем
            continue;
        if(!do_merge)
            t = new HandleTree(data[i], treeName, E);
        else{
            cout << "Merge, " << GetEnergyFromString(data[i]) << endl;
            cout << "WAS: " << t->getEntriesWithConditions() << endl; 
            t->Merge(data[i], GetEnergyFromString(data[i]));
            cout << "IS: " << t->getEntriesWithConditions() << endl;}
            
        t->makeFit();
        gPad->WaitPrimitive();
        cout << "Merge it?: ";
        cin >> do_merge;
        if(!do_merge||(i==n-1)){
            mdata.insert(mdata.end(), t);
            fits.insert(fits.end(), t->getFit());
            hists.insert(hists.end(), t->getHist());
        }
    }
    return;
}

double Dataset::GetLuminosity(double E){
    int n = lum.size();
    for(int i=0; i<n; i++){
        if( fabs(lum[i].first - E) < 0.00001 )
            return lum[i].second[0];
    }
    return 0;
}

vector<TF1*> Dataset::GetFits(){
    return fits;
}

vector<TH1D*> Dataset::GetHists(){
    return hists;
}

vector<TGraphAsymmErrors*> Dataset::GetTriggers(){
    int n = mdata.size();
    vector<double> E(n);
    vector<vector<double>> T(3);
    vector<vector<double>> C(3);
    vector<vector<double>> TC(3);
    double** res;
    for(int i=0; i<n; i++){
        res = mdata[i]->triggerEfficiency();
        E[i] = mdata[i]->getEnergy();
        T[0][i] = res[0][0];
        T[1][i] = res[1][0];
        T[2][i] = res[2][0];
        C[0][i] = res[0][1];
        C[1][i] = res[1][1];
        C[2][i] = res[2][1];
        TC[0][i] = res[0][2];
        TC[1][i] = res[1][2];
        TC[2][i] = res[2][2];
    }
    TGraphAsymmErrors* t = new TGraphAsymmErrors(n, &E[0], &T[0][0], 0, 0, &T[2][0], &T[1][0]);
    TGraphAsymmErrors* c = new TGraphAsymmErrors(n, &E[0], &C[0][0], 0, 0, &C[2][0], &C[1][0]);
    TGraphAsymmErrors* tc = new TGraphAsymmErrors(n, &E[0], &TC[0][0], 0, 0, &TC[2][0], &TC[1][0]);
    vector<TGraphAsymmErrors*> vg = {t, c, tc};
    return vg;
}

double Dataset::GetRadcor(double E){
    int n = radcor.size();
    double index[2] = {0, 1};
    double res;
    double x1, x2;
    double y1;
    double y2;
    
    for(int i=0; i<n-1; i++){
        
        if( (E>=radcor[i].first)&&(E<=radcor[i+1].first) ){
            index[0] = i;
            index[1] = i+1;
            break;
        }
    }
    if( E>=(radcor[n-1].first) ){
        index[0] = n-2; index[1] = n-1;
    }
    
    x1 = radcor[index[0]].first;
    x2 = radcor[index[1]].first;
    y1 = radcor[index[0]].second;
    y2 = radcor[index[1]].second;
    res = LinearApprox(E, x1, y1, x2, y2);
    return res;
}


TGraphAsymmErrors* Dataset::GetCS(){
    //gROOT->SetBatch(kTRUE);
    int n = mdata.size();
    vector<double> E(n);
    vector<vector<double>> cs(3);
    double** T;
    double N;
    double* R;
    double r;
    double L;
    cout << "DDD: " << n << endl;
    for(int i=0; i<n; i++){
        E[i] = mdata[i]->getEnergy();
        if(( mdata[i]->getEnergy() )<0 ){ //отрицательная энергия указывает на некорректно определённые светимости в чейне
            cout << "Warning!!!" << endl;
            cs[0][i] = 0;   cs[1][i] = 0;   cs[2][i] = 0;
            continue;}
        T = mdata[i]->triggerEfficiency();
        N = fits[i]->GetParameter(0);
        R = RegistrationEff(E[i]);
        r = GetRadcor(E[i]);
        L = mdata[i]->getLum();
        cout << N << '\t' << L << '\t' << T[0][2] << '\t' << r << '\t' << R[0] << '\t' << N/(L*T[0][2]*r*R[0]) << endl;
        cs[0].insert( cs[0].end(), N/(L*T[0][2]*r*R[0]) );
        cs[1].insert( cs[1].end(), cs[0][i]*sqrt( pow((fits[i]->GetParError(0))/N,2) + pow(T[1][2]/T[0][2],2) + pow(R[1]/R[0],2) ) );
        cs[2].insert( cs[2].end(), cs[0][i]*sqrt( pow((fits[i]->GetParError(0))/N,2) + pow(T[2][2]/T[0][2],2) + pow(R[2]/R[0],2) ) );
    }
    //gROOT->SetBatch(kFALSE);
    cross_sections = cs;
    TGraphAsymmErrors* t = new TGraphAsymmErrors(n, &E[0], &cs[0][0], 0, 0, &cs[2][0], &cs[1][0]);
    t->Draw();
    return t;
}

void Dataset::Save(string path = "Outputs/result.root"){
    TFile* f = TFile::Open(path.c_str(), "recreate");
    TTree* t = new TTree("t", "Tree");
    TF1 fitF;
    TH1D histF;
    double cs[3];
    double* regEff = NULL;
    double E;
    
    t->Branch("fitF", "TF1", &fitF, 32000, 0);
    //cout << sizeof(TH1D) << endl;
    t->Branch("histF", "TH1D", &histF, 32000, 0);
    t->Branch("E", &E, "E/D");
    t->Branch("cs", &cs, "cs[3]/D");
    //t->Branch("regEff", regEff, "regEff[3]/D");
    
    int n = fits.size();
    for(int i=0; i<n; i++){
        fitF = *fits[i];
        histF = *hists[i];
        cs[0] = cross_sections[0][i];
        cs[1] = cross_sections[1][i];
        cs[2] = cross_sections[2][i];
        E = model[i]->getEnergy();
        //regEff = RegistrationEff(E);
        t->Fill();
        //delete regEff;
    }
    
    t->Write();
    f->Close();
    //delete t;
    return;
}

void Dataset::Read(string path = "Outputs/result.root"){
    TFile* f = TFile::Open(path.c_str());
    fits.clear();
    hists.clear();
    cross_sections.clear();
    
    TTree* t = (TTree*)(f->Get("t"));
    
    TTreeReader theReader(t);
    TTreeReaderValue<TF1> fitBranch(theReader, "fitF");
    TTreeReaderValue<TH1D> histBranch(theReader, "histF");
    TTreeReaderArray<double> csBranch(theReader, "cs");

    while(theReader.Next()){
       TF1 f = *fitBranch;
       TH1D h = *histBranch;
       fits.insert(fits.end(), &f);
       hists.insert(hists.end(), &h);
       vector<double> cs0(3);
       cs0[0] = csBranch[0]; cs0[1] = csBranch[1]; cs0[2] = csBranch[2];
       cross_sections.insert(cross_sections.end(), cs0);
    }
    
    f->Close();
    return;
}

Dataset::~Dataset(){
    mdata.clear();
    fits.clear();
    hists.clear();
    data.clear();
    model.clear();
    lum.clear();
    radcor.clear();
    cross_sections.clear();
}
