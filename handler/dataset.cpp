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
    return ( (y2 - y1)/(x2 - x1) ) * (E - x1) + y1;
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
    gROOT->SetBatch(kTRUE);
    int n = model.size();
    double index[2] = {0, 1};
    double* res = new double [3];
    double x1, x2;
    double* y1;
    double* y2;
    
    for(int i=0; i<n-1; i++){
        
        if( (E>(model[i]->getEnergy()))&&(E<(model[i+1]->getEnergy())) ){
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
    gROOT->SetBatch(kFALSE);
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
    bool do_merge = false;
    for(int i=startPosition; i<n; i++){
        if(!do_merge)
            t = new HandleTree(data[i], treeName, GetEnergyFromString(data[i]));
        else
            t->Merge(data[i], GetEnergyFromString(data[i]));
            
        t->makeFit();
        cout << "Merge it?" << endl;
        cin >> do_merge;
        if((!do_merge)||(i==n-1)){
            mdata.insert(mdata.end(), t);
            fits.insert(fits.end(), t->getFit());
            hists.insert(hists.end(), t->getHist());
        }
        gPad->WaitPrimitive();
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
    for(int i=0; i<n; i++){
        if( fabs(radcor[i].first - E) < 0.00001 )
            return radcor[i].second;
    }
    cout << Form( "Radcor for %.2E isn't found", E) << endl; 
    return 1;
}

TGraphAsymmErrors* Dataset::GetCS(){ //дописать!!
    int n = mdata.size();
    vector<double> E;
    vector<vector<double>> cs(3);
    double** T;
    double N;
    double* R;
    double r;
    for(int i=0; i<n; i++){
        E[i] = mdata[i]->getEnergy();
        if(GetLuminosity(E[i])==0){
            cs[0][i] = 0;   cs[1][i] = 0;   cs[2][i] = 0;
            continue;}
        T = mdata[i]->triggerEfficiency();
        N = fits[i]->GetParameter(0);
        R = RegistrationEff(E[i]);
        r = GetRadcor(E[i]);
        cs[0][i] = N/(GetLuminosity(E[i])*T[0][2]*r*R[0]);
        cs[1][i] = cs[0][i]*sqrt( pow((fits[i]->GetParError(0))/N,2) + pow(T[1][2]/T[0][2],2) + pow(R[1]/R[0],2) );
        cs[2][i] = cs[0][i]*sqrt( pow((fits[i]->GetParError(0))/N,2) + pow(T[2][2]/T[0][2],2) + pow(R[2]/R[0],2) );
    }
    
    cross_sections = cs;
    TGraphAsymmErrors* t = new TGraphAsymmErrors(n, &E[0], &cs[0][0], 0, 0, &cs[2][0], &cs[1][0]);
    return t;
}

void Dataset::Save(string path = "Outputs/result.root"){
    TFile* f = TFile::Open(path.c_str(), "recreate");
    TTree* t = new TTree("t", "Tree");
    TF1* fitF = NULL;
    TH1D* histF = NULL;
    double cs[3];
    double* regEff = NULL;
    double E;
    t->Branch("fitF", "TF1", fitF, sizeof(TF1));
    t->Branch("histF", "TH1D", histF, sizeof(TH1D));
    t->Branch("cs", &cs, "cs[3]/D");
    t->Branch("regEff", regEff, "regEff[3]/D");
    t->Branch("E", &E, "E/D");
    
    int n = mdata.size();
    for(int i=0; i<n; i++){
        fitF = fits[i];
        histF = hists[i];
        cs[0] = cross_sections[i][0];
        cs[1] = cross_sections[i][1];
        cs[2] = cross_sections[i][2];
        E = model[i]->getEnergy();
        regEff = RegistrationEff(E);
        t->Fill();
        delete regEff;
    }
    t->Write();
    f->Close();
    delete t;
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
