#ifndef FIT_H
#define FIT_H

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TMath.h>
#include <TStyle.h>
#include <TTree.h>
#include <TChain.h>
#include <TGraphAsymmErrors.h>
#include <fstream>
#include <iostream>
#include <TTreeReader.h>
#include <string>
#include <vector>

const int n = 70;

using namespace TMath;
using namespace std;

class HandleTree{

    TChain* chain = NULL;
    TH1D* hist = NULL; //распределение по инвариантной массе
    TF1* fitF = NULL; //фит по инвариантной массе 
    std::vector<double> energies;
    std::vector<std::vector<double>> lums;
    string conditions = "win&m>450&m<550"; //условия на hist, fitF и всё остальное тоже
    double EL[3]; //сумма произведений энергий на светимости
    double L[3]; //полная светимость
    
    double variance(int, int);

public:

    HandleTree(string treeName){
        chain = new TChain(treeName.c_str());
    }
    
    HandleTree(string path, string treeName, std::vector<double> lum){   //path - путь до файла с деревом, treeName - имя дерева в рут файле
        std::cout << path.substr( path.find_last_of('/') + 2, path.length() - 5) << std::endl;
        //energies.insert( energies.end(), atof(path.substr( path.find_last_of('/') + 2, path.length() - 5).c_str()) );
        std::vector<double> zeros = {0,0,0};
        lums.insert( lums.end(), lum[0]>0 ? lum : zeros  );
        chain = new TChain(treeName.c_str());
        chain->Add(path.c_str());
        energies.insert( energies.end(), chain->GetMaximum("be") ); //be одинаковые, поэтому можно брать хоть min, хоть max
    }
    
    int* getTriggers(); //необходимо, чтоб в дереве переменная триггера называлась t и имела значение 0 - TF, 1 - CF, 2 - TF&CF
    vector<int> GetTriggers(); //новое поколение той же функции
    void Merge(string path, std::vector<double> lum); //путь к файлу другого дерева, дерево в этом файле должно называться так же как и в оригинальном
    double* getRegistrationEfficiency();
    double** triggerEfficiency();
    
    void makeHist(); //гистограмма по инвариантной массе ("m") для событий прошедших this->conditions
    void makeFit(); //сделать фит для гистограммы из getHist
    TH1D* getHist(){
        return this->hist;  }
    TF1* getFit(){
        return this->fitF;   }
    std::vector<double> getEnergy(); // с ошибками [E, +dE, -dE]
    std::vector<double> getLum(); // с ошибками [L, +dL, -dL]
    int getEntriesWithConditions(){
        return chain->GetEntries(conditions.c_str());   }
    
    
    ~HandleTree(){
        delete hist;
        delete chain;
        delete fitF;
        energies.clear();
        lums.clear();
        hist = NULL;
        chain = NULL;
        fitF = NULL;
    }

};

#endif
