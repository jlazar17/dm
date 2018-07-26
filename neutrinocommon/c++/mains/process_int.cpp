/*
Author  : C.A.Arguelles
Date    : 25/MAR/2013
*/

#include <vector>
#include <iostream>
#include "body.h"
#include "physconst.h"
#include "neuosc.h"
#include "neurho.h"
#include "interface.h"

#include "tools.h"
#include "xsections.h"

using namespace std;

int main(int argc, char **argv)
{
    string in_file_path;
    string out_file_path;

    if(argc != 3){
        printf("Usage : process.exe in_file_path out_file_path \n");
        exit(0);
    } else {
        in_file_path = argv[1];
        out_file_path = argv[2];
    }
    
    // setting up oscillation configuration
    // E_range and numneu will be setup from loaded file
    PhysConst* pc = new PhysConst();
    // setup oscillation
    OscSetup* setup = new OscSetup;
    setup->param = pc;
    // setup interactions
    setup->oscillation = true;
    setup->attenuation = true;
    setup->nc_inter = true;
    setup->cc_inter = false;
    setup->tau_regeneration = true;
    setup->muon = false;
    setup->electron = false;
    
    setup->neutype = 2;
    
    // put body and trayectory

    double phi_zenith = 3.14159265;
    EarthAtm* earth = new EarthAtm;
    setup->body = earth;
    setup->track = new EarthAtm::Track(phi_zenith);
    /*
    double baseline = 10000.0*pc->km;
    ConstantDensity* constdens = new ConstantDensity(20.0,0.5);
    ConstantDensity::Track* constdens_track = new ConstantDensity::Track(0.0,baseline);
    setup->body = constdens;
    setup->track = constdens_track;
    */
    VectorFluxState *initial_state_vector[2];
    
    InitLoadFlux(initial_state_vector,in_file_path,setup);
    
    VectorFluxState** final_state_vector = CalNeuOscInt(initial_state_vector,setup,1.0e-20,1.0e-20,true);
    //VectorFluxState** final_state_vector = CalNeuOscSUN(initial_state_vector,setup,1.0e-20,1.0e-20);

    // milk the results and print them.
    VectorFluxState* final_state_neu = final_state_vector[0];
    VectorFluxState* final_state_aneu = final_state_vector[1];
    
    ofstream out_file;
    out_file.open(out_file_path.c_str());
    
    for(unsigned int e = 0; e<setup->E_range.size();e++){
        out_file << setup->E_range[e] << "\t";
        
        for(int i = 0; i < setup->param->numneu; i ++){
            out_file << TraceFlavorProject(i,final_state_neu)->data[e] << "\t";            
        }
        for(int i = 0; i < setup->param->numneu; i ++){
            out_file << TraceFlavorProject(i,final_state_aneu)->data[e] << "\t";            
        }
        
        out_file << endl;
    }
    out_file.close();
    
    return 0;
}