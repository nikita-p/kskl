#ifndef DATASET_H
#define DATASET_H

#include "fit.h"
#include <vector>
#include <TSystem.h>
#include <TTreeReaderArray.h>
#include <TPad.h>
#include <fstream>
#include <utility>
#include <string>
#include <iostream>
#include <algorithm>


using namespace std;

class Dataset{

    vector<vector<double>> cross_sections;

    vector<HandleTree*> mdata;
    vector<TF1*> fits;
    vector<TH1D*> hists;
    
    string modelFolder;
    string dataFolder;
    string lumFile;
    string radcorFile;
    
    vector<string> data;
    vector<HandleTree*> model;
    vector<pair<double, double*>> lum;
    vector<pair<double, double>> radcor;
    
    string treeName = "InvMass";
    
    double* LinearApprox(int, double, double, double*, double, double*);
    double LinearApprox(double, double, double, double, double);
    double GetEnergyFromString(string);
    static bool compt(HandleTree*, HandleTree*);
    static bool comps(string, string);
    
    void GetDataVector();
    void GetModelVector();
    void GetLumVector();
    void GetRadcorVector();
    
    void ClearFits();
    void ClearHists();
    void ClearMdata();

public:
    Dataset(string model, string data, string lum, string radcor);
    double* RegistrationEff(double);
    double GetLuminosity(double);
    double GetRadcor(double);
    void AutoFit(); // автоматический пакетный фит гистограмм без мёрджа
    void PowerFit(int); //фитировать гистограммы по очереди (можно даже руками), можно мёрджить
    
    vector<TF1*> GetFits();
    vector<TH1D*> GetHists();
    vector<TGraphAsymmErrors*> GetTriggers();
    TGraphAsymmErrors* GetEffRegs(); //написать
    TGraphAsymmErrors* GetLums(); //написать
    TGraphAsymmErrors* GetCS();
    
    void Save(string);
    void Read(string);
    ~Dataset();
};


#endif
