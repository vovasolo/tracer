#include <iostream>
#include <fstream>
#include <vector>
#include "tracer.h"

void writeColumnData(vector<vector<float> > *data, std::ofstream &file);

int main(int argc, char argv[])
{
    vector<vector<float> > result;

    Tracer *t = new Tracer();
    t->SetLogging(false);
    
    
    std::cout << "Loading geometry" << std::endl;
    t->Init("geometry/geo.root", "geometry/ovr.txt");
    std::cout << "Overrides:" << std::endl;
    t->PrintOmap();
    t->SetMaxSteps(1000);
    
//    t->OnePhoton(0,0,0);
//    t->PhotonBomb(0,0,0,100);
    t->ReadScanData("data/scan.txt");
    t->Scan(1000, result);
    
    std::ofstream fileout("out.txt");
    writeColumnData(&result, fileout);
    fileout.close();

    std::cout << "Done! Results are in out.txt" << std::endl;
    return 0;
}
