#ifndef DATASET_H
#define DATASET_H

#include "fit.h"
#include <vector>
#include <TSystem.h>
#include <fstream>
#include <utility>
#include <string>
#include <iostream>


using namespace std;

class Dataset{

    vector<TF1*> fits;
    vector<TH1D*> hists;
    
    string modelFolder;
    string dataFolder;
    string lumFile;
    string radcorFile;
    
    vector<HandleTree*> data;
    vector<HandleTree*> model;
    vector<pair<double, double*>> lum;
    vector<pair<double, double>> radcor;
    
    string treeName = "InvMass";
    
    double* LinearApprox(int, double, double, double*, double, double*);
    double LinearApprox(double, double, double, double, double);
    static bool comp(HandleTree*, HandleTree*);
    
    void GetDataVector();
    void GetModelVector();
    void GetLumVector();
    void GetRadcorVector();

public:
    Dataset(string, string, string, string);
    double* RegistrationEff(double);
    void AutoFit();
    
    vector<TF1*> GetFits();
    vector<TH1D*> GetHists();
};


#endif
