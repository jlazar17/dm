#ifndef __NEURHO_H
#define __NEURHO_H

#include "global.h"
#include "xsections.h"
#include "taudecay.h"
#include "body.h"
#include "physconst.h"
#include "neuosc.h"
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
//#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_odeiv2.h>


//#define AttenuationOnlyCC

//#define InitInter_DEBUG
#define InitInter_Unitarity_FIX
//#define InitInter_Unitarity_DEBUG
//#define UpdateInter_DEBUG
//#define CalNeuOscInt_DEBUG
//#define CalNeuOscIntStep_DEBUG
//#define ConvertVectorStateToVector_DEBUG
//#define RHS_RHO_DEBUG
//#define RHS_Hamiltonian_DEBUG
//#define RHS_Matrix_DEBUG
//#define RHS_Matrix_Osc_DEBUG
//#define RHS_Matrix_Int_DEBUG
//#define RHS_Matrix_NC_DEBUG
//#define RHS_Matrix_IntVector_DEBUG
//#define RHS_Matrix_Tau_DEBUG
//#define RHS_Matrix_tau_lepton_DEBUG
//#define RHS_Matrix_tau_neutrino_DEBUG

#define RK_Unitarity_DEBUG

//#define RHS_Matrix_Unitarity_DEBUG

//#define CUDA_ENABLE
// CUDA vector addition is slower
// CPU vector addition

//#define USE_GSL_ODE_DRIVER

#define ManualTauReinjection
//#define TauRegManual_DEBUG
//#define TauLeptonReinjection_DEBUG

//#define USE_FIX_STEP


typedef struct
{
    // energy values
    public :
    vector<double> E_range;
    Track* track;
    PhysConst* param;

    Body* body;
    // control interaction types
    bool oscillation;
    bool attenuation;
    bool nc_inter;
    bool cc_inter;
    bool tau_regeneration;
    bool muon;
    bool electron;
    
    int neutype;
    int basis;
    // mixing matrix and useful to know matrices
    gsl_matrix_complex* flavorM2;
    gsl_matrix_complex* massM2;
    gsl_matrix_complex* U;
    // vector proyector
    vector<gsl_matrix_complex*> FlavorProjector;
    vector< vector<gsl_matrix_complex*> > VectorFlavorProjector;
    // useful memory matrices
    gsl_matrix_complex* RHS_matrix;
    
} OscSetup;

typedef struct
{
    // This structure contains all the
    // information to build the RHS
    // interaction terms
    public :    
    // decay lengths
    vector<gsl_vector_complex*> Lint;
    vector<gsl_vector_complex*> LNC;
    vector<gsl_vector_complex*> LCC;
    gsl_vector_complex* LTauDecay;

    // neutrino cross sections
    vector<gsl_vector_complex*> sigma_CC;
    vector<gsl_vector_complex*> sigma_NC;
    
    // neutrino dif. cross sections
    vector<gsl_matrix_complex*> dsignudE_NC;
    vector<gsl_matrix_complex*> dsignudE_CC;
    
    // lepton decay dif. cross sections
    vector<gsl_matrix_complex*> dNledE_NC;
    vector<gsl_matrix_complex*> dNledE_CC;
    
    // tau decay spectra
    gsl_matrix_complex* dNdENuTau_All;
    gsl_matrix_complex* dNdENuTau_Lep;
    gsl_matrix_complex* dNdELeTau_All;
    gsl_matrix_complex* dNdELeTau_Lep;
} InterOperator;

typedef struct
{
    public :
    OscSetup* setup;
    InterOperator* interactions;
    int numeqn;
} Box;

typedef struct
{
    public :
    // fluxes
    gsl_matrix_complex* F_nu;
    gsl_matrix_complex* F_le;
    double* FluxFlavor;
    double* FluxMass;
    SU_vector FluxState;
} FluxState;

typedef vector<FluxState*> VectorFluxState;

// tools

gsl_matrix_complex* Projector(int,int);
gsl_matrix_complex* FlavorProjector(int,int);
gsl_matrix_complex* FlavorProjector(int,OscSetup*);
vector<gsl_matrix_complex*> InitProjectorVector(OscSetup*);
gsl_vector* TraceFlavorProject(int,VectorFluxState*);

InterOperator* InitInter(OscSetup*);
int UpdateInter(InterOperator*,OscSetup*);

// flux initializer
VectorFluxState* InitSimpleState(int,OscSetup*);
VectorFluxState* InitPowerLawFluxState(int,double,OscSetup*);
VectorFluxState* InitBodyFluxState(OscSetup*);

// vector state converters
gsl_vector_complex* ConvertVectorStateToVector(VectorFluxState*,OscSetup*);

VectorFluxState* ConvertVectorToVectorState(gsl_vector_complex*,OscSetup*);

// tau reinjection algorithms
int TauLeptonReinjection(double *, void *);

int TauLeptonReinjection(double *, void *, double *, void *);

// calculation
gsl_matrix_complex* CalculateRHSMatrix(double,InterOperator*,OscSetup*);

int RHS_RHO(double,const double *, double * ,void *);
int JAC_RHO(double,const double *, double * , double [], void *);

// runge kutta implementation
void ConfigureSetup(OscSetup*);
VectorFluxState** CalNeuOscInt(VectorFluxState**,OscSetup*,double,double,bool);

#endif