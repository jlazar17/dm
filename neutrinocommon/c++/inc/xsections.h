#ifndef __XSECTIONS_H
#define __XSECTIONS_H

#include "global.h"
#include "tools.h"
#include "physconst.h"
#include <string>
#include <cmath>
#include <math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

//#define NCS_Init_DEBUG

/*
This functions return the neutrino nucleon
assuming isoscalar target charged (CC)
and neutral (NC) deep inelastic
cross sections obtained from NUSIGMA.
*/

// NUSIGMA FORTRAN ROUTINES
extern "C" double dsdxdy_(double*,double*,double*,int*,int*,int*);
extern "C" double dsde_(double*,double*,int*,int*,int*);
extern "C" double sigma_(double*,int*,int*,int*);

// neutrino cross sections

class NeutrinoCrossSections{
    public :
        
    double Emin;
    double Emax;
    int div;
        
    Table dsde_CC_tbl;
    Table dsde_NC_tbl;
    Table sigma_CC_tbl;
    Table sigma_NC_tbl;
    
    NeutrinoCrossSections(double,double,int);
    
    double dsde_CC(int,int,int,int);
    double dsde_NC(int,int,int,int);
    double sigma_CC(int,int,int);
    double sigma_NC(int,int,int);
    
    /*
    double Calc_dsde_CC(double,double,int,int);
    double Calc_dsde_NC(double,double,int,int);
    double Calc_sigma_CC(double,int,int);
    double Calc_sigma_NC(double,int,int);
    */
};

// tau lepton cross sections


#endif
