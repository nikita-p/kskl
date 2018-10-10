#ifndef COPIER_H
#define COPIER_H

#include "fit.h"
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TSystem.h>

#include <vector> //for vector
#include <iterator> //for vector
#include <string> //for string
#include <algorithm>

using namespace std;

class Copier{
  
  string folder;
  vector<HandleTree*> vec;
  
  
  static bool comp (HandleTree* a, HandleTree* b) {
    return (a->getEnergy()) < (b->getEnergy());
  }
    
  public:
    
    Copier(string folder){ //абсолютный путь до папки
        void* dir = gSystem->OpenDirectory( folder.c_str() );
    
        if( !dir ){
            std::cout << "Bad times" << std::endl;
            return; }
        
        const char* entry;
        string end = ".root";
        TString filename;
        
        entry = gSystem->GetDirEntry(dir);
        
        while( entry ){
            filename = entry;
            if(filename.EndsWith(end.c_str())){
                filename = gSystem->ConcatFileName(folder.c_str(), entry);
                std::cout << filename << std::endl;
                //загружать файл и в функцию
                HandleTree* f = new HandleTree(string(filename),"InvMass");
                vec.insert(vec.end(), f);
            }
            entry = gSystem->GetDirEntry(dir);
        }
        sort(vec.begin(), vec.end(), comp);
    }

    vector<TH1D*> showHists();
    vector<TF1*> showFits();
    vector<TGraphAsymmErrors*> showTriggers();
    void writeInformation(string);

};


#endif
