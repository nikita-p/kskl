#include "Trph.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

#include <iostream>
#include <string>
#include <ctime>
#include <fstream>
#include <TTree.h>
#include <TFile.h>

#define mKs 497.614
#define mPi 139.570

using namespace std;

class TreeReader{
public:

    TreeReader(){
    }
    
    void washing(string inputPath, string outputPath, bool time=false){
        double time0 = clock()/(CLOCKS_PER_SEC+0.);
        TFile *input = TFile::Open(inputPath.c_str());
        TFile *output = TFile::Open(outputPath.c_str(), "recreate");
        TTree *tr = (TTree*)(input->Get("tr_ph"));
        
        Trph cl(tr);
        cl.Loop();
        
        output->Write();
        output->Close();
        
        if(time) cout << "Time for washing: " << clock()/(CLOCKS_PER_SEC+0.) - time0 << "seconds." << endl;
        delete output;
        return;    
    }
    
    void washingFromFile(string filePath){
        ifstream f(filePath.c_str());
        string input, output;
        while( !f.eof() ){
            f >> input >> output;
            if(input=="")
                break;
            washing(input, output, true);
            cout << "Well done" << endl << input << '\t' << output << endl;
        }
        return;
    }
    
};


void events(string year){
    TreeReader t;
    if(year=="11"||year=="2011"){
        cout << "Go, " << year << endl;
        t.washingFromFile("Inputs/trees11.dat");    }
    if(year=="12"||year=="2012"){
        cout << "Go, " << year << endl;
        t.washingFromFile("Inputs/trees12.dat");    }
    if(year=="17"||year=="2017"){
        cout << "Go, " << year << endl;
        t.washingFromFile("Inputs/trees17.dat");    }
    if(year=="model"){
        cout << "Go, " << year << endl;
        t.washingFromFile("Inputs/treesModel.dat");    }
    return;
}
