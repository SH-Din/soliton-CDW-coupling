#include "parameters.hpp"
#include "free_energy_construction.hpp"
#include <fstream>

int main(){
    try{
        // try to read parameters.json; fall back to defaults if missing
        Parameters pars;
        std::ifstream parsFile("parameters.json");
        if(parsFile.good()) pars = Parameters("parameters.json");
        pars.print();
        FreeEnergyConstruction fec(pars);
        fec.cal();
    }catch(const std::exception& e){
        std::cerr<<"Error: "<<e.what()<<"\n";
        return 1;
    }
    return 0;
}