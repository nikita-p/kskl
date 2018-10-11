#include "dataset.h"

using namespace std;

bool Dataset::comp(HandleTree* a, HandleTree* b){
    return (a->getEnergy())<(b->getEnergy());
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
            HandleTree* f = new HandleTree(string(filename), treeName.c_str());
            data.insert(data.end(), f);
        }
        entry = gSystem->GetDirEntry(dir);
    }
    sort(data.begin(), data.end(), comp);
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
            HandleTree* f = new HandleTree(string(filename), treeName.c_str());
            model.insert(model.end(), f);
        }
        entry = gSystem->GetDirEntry(dir);
    }
    sort(model.begin(), model.end(), comp);
    return;
}

void Dataset::GetLumVector(){
    ofstream flum(lumFile.c_str());
    
    double e;
    while(!flum.eof()){
        double* l = new double [3];
        flum << e << l[0] << l[1];
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
    ofstream frad(radcorFile.c_str());
    
    double e;
    double rad;
    while(!frad.eof()){
        frad << e << rad;
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
    return res;
}

void Dataset::AutoFit(){
    
    int n = data.size();
    for(int i=0; i<n; i++){
        data[i]->makeFit();
        fits.insert(fits.end(), data[i]->getFit());
        hists.insert(hists.end(), data[i]->getHist());
    }
    return;
}

vector<TF1*> Dataset::GetFits(){
    return fits;
}

vector<TH1D*> Dataset::GetHists(){
    return hists;
}


