#include "dataset.h"
#include <fstream>

void SaveToFile(Dataset* a){
    ofstream o("Outputs/result.txt");
    o.write((char*)a, sizeof(Dataset));
    o.close();    
    return;
}

void ReadFromFile(Dataset* a){
    ifstream o("Outputs/result.txt");
    o.read((char*)a, sizeof(Dataset));
    o.close();    
    return;
}
