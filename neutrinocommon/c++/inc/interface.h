#ifndef __INTERFACE_H
#define __INTERFACE_H

#include "global.h"
#include "body.h"
#include "physconst.h"
#include "neuosc.h"
#include "neurho.h"
#include "neualg.h"
#include "tools.h"
//#define Mask_DEBUG
//#define CalNeuOscIntMask_DEBUG
#define CalNeuOscSUNMask_DEBUG

struct BodyTrack {
    Body* body;
    Track* track;
};

struct FluxBox{
    double* E_values;
    
    double* nue_values;
    double* numu_values;
    double* nutau_values;
    
    double* anue_values;
    double* anumu_values;
    double* anutau_values;
    
    FluxBox(int size){
        E_values = new double[size];
        
        nue_values = new double[size];
        numu_values = new double[size];
        nutau_values = new double[size];
        
        anue_values = new double[size];
        anumu_values = new double[size];
        anutau_values = new double[size];
    }
};

// flux interface
int InitLoadFlux(VectorFluxState**,string,OscSetup*);
int InitLoadFlux(VectorFluxState**,Table,OscSetup*);
VectorFluxState** InitLoadFlux(string,OscSetup*);
VectorFluxState* InitGeneralFlux(FluxBox*,OscSetup*,int);


// RK interface
int SetOscParamMask(double [],PhysConst*);
int SetBodyTrackMask(int, double [],double [],Body*,Track*);
BodyTrack* GetBodyTrackMask(int, double [],double []);

vector<double> CalNeuOscGSLMask(int,double,int,double [],double [],int,double [],double,double,int,int);

vector< vector<double> > CalNeuOscIntMask(int, double [],
                                        double,double,int,
                                        int,double [],double [],
                                        int,double [],
                                        double,double,
                                        int,int,int,int,
                                        int,int,int,
                                        int);

vector< vector<double> > CalNeuOscSUNMask(int, double [],
                                        double,double,int,
                                        int,double [],double [],
                                        int,double [],
                                        double,double,
                                        int,int,int,int,
                                        int,int,int,
                                        int);

#endif