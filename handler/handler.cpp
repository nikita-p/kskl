#include "fit.h"
#include <TSystem.h>
#include <TString.h>

void handler(string mode){
    
    mode = "Outputs/trees/"+mode;
    void* dir = gSystem->OpenDirectory( mode.c_str() );
    
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
            filename = gSystem->ConcatFileName(mode.c_str(), entry);
            std::cout << filename << std::endl;
            //загружать файл и в функцию
            HandleTree f(string(filename),"InvMass");
            std::cout << f.getRegistrationEfficiency()[0] << "\t" << f.getRegistrationEfficiency()[1] << endl;        
        }
        
        entry = gSystem->GetDirEntry(dir);
    }
    
    HandleTree *f = new HandleTree("Outputs/trees/2011/e550.root", "InvMass");
    f->makeFit();
    
    f->getHist()->Draw();
    f->getFit()->Draw("same");
    
    return;
}
