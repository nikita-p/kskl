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

const int n = 70;

using namespace TMath;


class HandleTree{

    TChain* chain = NULL;
    TH1D* hist = NULL; //распределение по инвариантной массе
    TF1* fitF = NULL; //фит по инвариантной массе 
    double energy; //!!!поправить, т.к. не будет работать при мёрдже
    string conditions = "win&m>450&m<550"; //условия на hist, fitF и всё остальное тоже
    
    double variance(int, int);

public:

    HandleTree(string treeName){
        chain = new TChain(treeName.c_str());
    }
    
    HandleTree(string path, string treeName){   //path - путь до файла с деревом, treeName - имя дерева в рут файле
        std::cout << path.substr( path.find_last_of('/') + 2, path.length() - 5) << std::endl;
        energy = atof( path.substr( path.find_last_of('/') + 2, path.length() - 5).c_str() );
        chain = new TChain(treeName.c_str());
        chain->Add(path.c_str());
    }
    
    int* getTriggers(); //необходимо, чтоб в дереве переменная триггера называлась t и имела значение 0 - TF, 1 - CF, 2 - TF&CF
    void Merge(string path); //путь к файлу другого дерева, дерево в этом файле должно называться так же как и в оригинальном
    double* getRegistrationEfficiency();
    double** triggerEfficiency();
    
    void makeHist(); //гистограмма по инвариантной массе ("m") для событий прошедших this->conditions
    void makeFit(); //сделать фит для гистограммы из getHist
    TH1D* getHist(){
        return this->hist;  }
    TF1* getFit(){
        return this->fitF;   }
    double getEnergy(){
        return this->energy;    }
    
    
    ~HandleTree(){
        delete hist;
        delete chain;
        delete fitF;
        hist = NULL;
        chain = NULL;
        fitF = NULL;
    }

};

#endif
