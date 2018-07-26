#ifndef __NEUALG_H
#define __NEUALG_H

#include "global.h"
#include "xsections.h"
#include "taudecay.h"
#include "body.h"
#include "physconst.h"
#include "neuosc.h"
#include "neurho.h"
#include "SUNalg.h"

#ifdef GPU_ENABLE
    #include "gpu_tools.h"
#endif

#include <iostream>
#include <float.h>
#include <math.h>
#include <complex>
#include <vector>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_odeiv2.h>

#define ESIZEMAX 351
#define NEUINT 3 // number of active neutrino flavors
#define NEUMAX 4 // maximum number of neutrinos

#define USE_GSL_ODE_DRIVER
//#define CalNeuOscSUN_DEBUG
//#define CalNeuOscSUNStep_DEBUG
//#define UpdateInteractions_DEBUG
//#define RHS_SUN_DEBUG
//#define RHS_SUN_POS_DEBUG
//#define RHS_SUN_OSC_DEBUG
//#define RHS_SUN_ATT_DEBUG
//#define RHS_SUN_OUT_DEBUG
//#define RHS_SUN_NC_DEBUG
//#define RHS_SUN_TAU_DEBUG
//#define RHS_SUN_NC_EXT_DEBUG
#define RHS_SUN_OLD_ATT

//#define TauLeptonReinj_SUN_DEBUG

#define ManualTauReinjection
#define FixCrossSections
#define SterileNC

// master structures
typedef struct
{
    public :
    // basic
    int numeqn;
    double E_range[ESIZEMAX];
    double delE[ESIZEMAX];
    double br_lep;
    int esize;
    int neutype;
    int numneu;
    Track* track;
    Body* body;
    PhysConst* param;
    // control interaction types
    bool oscillation;
    bool attenuation;
    bool nc_inter;
    bool cc_inter;
    bool tau_regeneration;
    // cross section and lenghts
    double sigma_CC[NEUINT][ESIZEMAX];
    double sigma_NC[NEUINT][ESIZEMAX];
    double dNdE_CC[NEUINT][ESIZEMAX][ESIZEMAX];
    double dNdE_NC[NEUINT][ESIZEMAX][ESIZEMAX];
    double invlen_int[NEUINT][ESIZEMAX];
    double invlen_CC[NEUINT][ESIZEMAX];
    double invlen_NC[NEUINT][ESIZEMAX];
    double invlen_tau[ESIZEMAX];
    double dNdE_tau_all[ESIZEMAX][ESIZEMAX];
    double dNdE_tau_lep[ESIZEMAX][ESIZEMAX];
    // su_vectors
    SU_vector mass_projectors[NEUMAX];
    SU_vector flavor_projectors[NEUMAX];
    SU_vector evol_mass_proj[NEUMAX*ESIZEMAX];
    SU_vector evol_flavor_proj[NEUMAX*ESIZEMAX];
    SU_vector DM2;
} SUBox;

typedef struct
{
    public :
    // fluxes
    //vector<SU_vector> neu_array;
    //vector<double> lep_array;
    SU_vector neu_array[ESIZEMAX];
    double lep_array[ESIZEMAX];   
} SUState;

typedef struct
{
    public :
    vector<SU_vector> su_array;
    vector<double> flavor_array;
    vector<double> mass_array;
} NeutrinoFlux;

void InitSUBox(SUBox*);
void InitSUBoxNoInt(SUBox*);
void InitProjectors(SUBox*);
void InitHamiltonian(SUBox*);
void InitInteractions(SUBox*);
void UpdateInteractions(SUBox*);

//void RotateToMass(SU_vector*,PhysConst*);
//void RotateToFlavor(SU_vector*,PhysConst*);
void EvolveProjectors(double,SUBox*);

SU_vector SUHInt(int,SUBox*);

void TauLeptonReinj_SUN(double*,void*);
void TauLeptonReinj_SUN(double*,void*,double*,void*);
int RHS_SUN(double,const double *, double * ,void *);
int JAC_SUN(double,const double *, double * , double [], void *);

// some convertors 
void ConvertOscSetupToSUBox(OscSetup*,SUBox*);

SUState* ConvertVectorStateToSUState(VectorFluxState*,SUBox*);
SUState* ConvertSUVectorToSUState(SUVector*,SUBox*);
VectorFluxState* ConvertSUStateToVectorState(SUState*,SUBox*);

//VectorFluxState** CalNeuOscSUN(VectorFluxState**,OscSetup*,double,double);
int CalNeuOscSUN(VectorFluxState**,OscSetup*,double,double);
int CalNeuOscSUN(VectorFluxState**,SUBox**,double,double);
int CalNeuOscSUN(SUVector**,SUBox**,double,double);
#endif
