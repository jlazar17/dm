#ifndef __SUNALG_H
#define __SUNALG_H

#include <vector>
#include <cstring>
#include <iostream>
#include <float.h>
#include <math.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>

#include "global.h"
#include "physconst.h"

#define NEUMAXSQ 16 // related to maximum number of neutrinos

#define SQR(x)      ((x)*(x))

class SU_vector{
  public:
    int dim;
    //double* components;
    double components[NEUMAXSQ];
    
    // declaracion 
    SU_vector();
    SU_vector(int);
    SU_vector(int,double*);
    SU_vector(std::vector<double>);
    SU_vector(gsl_matrix_complex*);
    // initialization
    void InitSU_vector(gsl_matrix_complex*);
    void InitSU_vector(int);
    // destruction
    ~SU_vector();
    
    // operations
    SU_vector Rescale(double);
    SU_vector Rotate(int,int,double, double);
    SU_vector RotateToMass(const PhysConst*);

    bool operator ==(const SU_vector&);
    double operator*(const SU_vector&);
    
    SU_vector operator*(const double);
    SU_vector operator =(const SU_vector&);
    SU_vector operator+=(const SU_vector&);
    SU_vector operator-=(const SU_vector&);
    SU_vector operator +(const SU_vector&);
    SU_vector operator -(const SU_vector&);

    //ostream overload operator
    friend ostream& operator<<(ostream&, const SU_vector&);
};

SU_vector operator *(const double x, SU_vector&);

SU_vector SUEvolve(SU_vector,SU_vector,double);
SU_vector SUiConmutator(SU_vector,SU_vector);
SU_vector SUAnticonmutator(const SU_vector,const SU_vector);

double SUTrace(SU_vector*,SU_vector*);
double SUTrace(SU_vector,SU_vector);

void SUPrint(const SU_vector*);
void SUPrint(const SU_vector);

typedef vector< vector < SU_vector > > SUVectorTable;
typedef vector<SU_vector> SUVector;

#endif
