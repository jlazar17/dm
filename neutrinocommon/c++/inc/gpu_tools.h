#ifndef __GPU_TOOLS_H
#define __GPU_TOOLS_H

#include "global.h"
#include <gsl/gsl_matrix.h>
#include <cuda_runtime.h>
#include <iostream>
#include "cublas_v2.h"

//#define CUDAComplexMatrixAdd_DEBUG

using namespace std;

int CUDAComplexMatrixAdd(gsl_matrix_complex* A,gsl_matrix_complex* B);

#endif

