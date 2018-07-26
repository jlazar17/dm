/*
Author  : C.A.Arguelles
Date    : 14/DEC/2012

Test functions. Sandbox file.

*/

#include <vector>
#include <iostream>
#include "body.h"
#include "physconst.h"
#include "neuosc.h"
#include "neurho.h"
#include "neualg.h"

#include "tools.h"
#include "xsections.h"

using namespace std;

int main()
{
    cout << "TESTING DENSITY MATRIX OSCILLATIONS" << endl;
    PhysConst* pc = new PhysConst();
    pc->numneu = 3;
    pc->neutype = 0;

    //pc->muon = 0;
    //pc->tau = 1;
    //pc->th12 = 3.14159265/4.0;
    //pc->dm21sq = 2.6e-3

    pc->th12 = 0.563942;
    pc->th13 = 0.154085;
    pc->th23 = 0.785398;

    pc->dm21sq = 7.65e-05;
    pc->dm31sq = 0.00247;
    pc->delta1 = 0.0;

    pc->Refresh();

    // setup oscillation
    OscSetup* setup = new OscSetup;
    setup->param = pc;
    //setup->E_range = logspace(1.0e2*pc->GeV,1.0e6*pc->GeV,149);
    setup->E_range = logspace(1.0*pc->GeV,1.0e1*pc->GeV,59);
    // setup interactions
    setup->oscillation = true;
    setup->attenuation = false;
    setup->nc_inter = false;
    setup->cc_inter = false;
    setup->tau_regeneration = false;
    setup->muon = false;
    setup->electron = false;

    setup->neutype = 0;
    setup->basis = 2;

    /*
    cout << "Vacuum" << endl;  
    setup->body = new Vacuum;
    setup->track = new Vacuum::Track(0.0,500.0*pc->km);
    */

    /*
    cout << "Earth" << endl;
    Earth* earth = new Earth;
    double L = 500.0*pc->km;
    setup->body = earth;
    setup->track = new Earth::Track(0,L,L);
    */
    /*
    cout << "Sun" << endl;
    Sun* sun = new Sun;

    setup->body = sun;
    setup->track = new Sun::Track(0.0,sun->radius*pc->km);
    */

    cout << "EarthAtm" << endl;
    EarthAtm* earth = new EarthAtm;
    setup->body = earth;
    //setup->track = new EarthAtm::Track(3.14159265359);
    setup->track = new EarthAtm::Track(acos(-1.0));

    //VectorFluxState* initial_state = InitSimpleState(pc->muon,setup);
    double power = 0.0;//-2.0;//-1.0;

    VectorFluxState* initial_state = InitPowerLawFluxState(pc->muon,power,setup);
    VectorFluxState* initial_state_cp = InitPowerLawFluxState(pc->muon,power,setup);

    VectorFluxState *initial_state_vector[2];
    initial_state_vector[0] = initial_state;
    initial_state_vector[1] = initial_state_cp;

    CalNeuOscSUN(initial_state_vector,setup,1.0e-20,1.0e-20);
    VectorFluxState** final_state_vector = initial_state_vector;
    //VectorFluxState** final_state_vector = CalNeuOscInt(initial_state_vector,setup,1.0e-20,1.0e-20);

    VectorFluxState* final_state = final_state_vector[0];
    //VectorFluxState* final_state_cp = final_state_vector[1];
    /*
    cout << "INITIAL STATE" << endl;
    for(unsigned int e = 0; e<setup->E_range.size();e++){
        cout << "p_e E = " << setup->E_range[e]/pc->GeV << " [GeV] : " << (TraceFlavorProject(pc->electron,initial_state_cp)->data[e])/pow(setup->E_range[e],-1) << endl;
        cout << "p_mu E = " << setup->E_range[e]/pc->GeV << " [GeV] : " << (TraceFlavorProject(pc->muon,initial_state_cp)->data[e])/pow(setup->E_range[e],-1) << endl;
        cout << "p_tau E = " << setup->E_range[e]/pc->GeV << " [GeV] : " << (TraceFlavorProject(pc->tau,initial_state_cp)->data[e])/pow(setup->E_range[e],-1) << endl;
    }

    cout << "FINAL STATE" << endl;
    for(unsigned int e = 0; e<setup->E_range.size();e++){
        cout << "p_e E = " << setup->E_range[e]/pc->GeV << " [GeV] : " << (TraceFlavorProject(pc->electron,final_state)->data[e])/pow(setup->E_range[e],-1) << endl;
        cout << "p_mu E = " << setup->E_range[e]/pc->GeV << " [GeV] : " << (TraceFlavorProject(pc->muon,final_state)->data[e])/pow(setup->E_range[e],-1) << endl;
        cout << "p_tau E = " << setup->E_range[e]/pc->GeV << " [GeV] : " << TraceFlavorProject(pc->tau,final_state)->data[e]/pow(setup->E_range[e],-1) << endl; 
    }
    */

    /*
    cout << "INITIAL STATE" << endl;
    for(unsigned int e = 0; e<setup->E_range.size();e++){
        cout << "p_e E = " << setup->E_range[e]/pc->GeV << " [GeV] : " << (TraceFlavorProject(pc->electron,initial_state_cp)->data[e]) << endl;
        cout << "p_mu E = " << setup->E_range[e]/pc->GeV << " [GeV] : " << (TraceFlavorProject(pc->muon,initial_state_cp)->data[e]) << endl;
        cout << "p_tau E = " << setup->E_range[e]/pc->GeV << " [GeV] : " << (TraceFlavorProject(pc->tau,initial_state_cp)->data[e]) << endl;
    }

    cout << "FINAL STATE" << endl;
    for(unsigned int e = 0; e<setup->E_range.size();e++){
        cout << "p_e E = " << setup->E_range[e]/pc->GeV << " [GeV] : " << (TraceFlavorProject(pc->electron,final_state)->data[e]) << endl;
        cout << "p_mu E = " << setup->E_range[e]/pc->GeV << " [GeV] : " << (TraceFlavorProject(pc->muon,final_state)->data[e]) << endl;
        cout << "p_tau E = " << setup->E_range[e]/pc->GeV << " [GeV] : " << TraceFlavorProject(pc->tau,final_state)->data[e] << endl; 
    }
    */

    for(unsigned int e = 0; e<setup->E_range.size();e++){
        cout << setup->E_range[e]/pc->GeV << " " << (TraceFlavorProject(pc->electron,final_state)->data[e])/(TraceFlavorProject(pc->muon,initial_state_cp)->data[e]) << endl;
    }

}
