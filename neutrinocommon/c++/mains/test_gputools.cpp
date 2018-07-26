#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "gpu_tools.h"
#include <iostream>

int main(){
    int size = 3;
    
    gsl_matrix_complex* A = gsl_matrix_complex_calloc(size,size);
    gsl_matrix_complex* B = gsl_matrix_complex_calloc(size,size);
    
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            gsl_matrix_complex_set(A,i,j,gsl_complex_rect(pow(double(i+1),2.5),pow(double(j+1),-2)));
            gsl_matrix_complex_set(B,i,j,gsl_complex_rect(pow(double(i+1),3.5),pow(double(j+1),-3)));
        }
    }
    /*
    double* A_dbl = A->data;
    for (int i = 0; i < 2*size*size; i++){
        cout << A_dbl[i] << endl;
    }
    
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            std::cout << gsl_matrix_complex_get(A,i,j).dat[0] << " " << gsl_matrix_complex_get(A,i,j).dat[1] << std::endl;
         }
    }
    */
    CUDAComplexMatrixAdd(A,B);

    gsl_matrix_complex* C = gsl_matrix_complex_calloc(size,size);
    gsl_matrix_complex* D = gsl_matrix_complex_calloc(size,size);
    
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            gsl_matrix_complex_set(C,i,j,gsl_complex_rect(pow(double(i+1),2.5),pow(double(j+1),-2)));
            gsl_matrix_complex_set(D,i,j,gsl_complex_rect(pow(double(i+1),3.5),pow(double(j+1),-3)));
        }
    }

    gsl_matrix_complex_add(C,D);
    
    /*
    gsl_matrix_complex_sub(A,C);
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            std::cout << gsl_matrix_complex_get(A,i,j).dat[0] <<
            "+i" << gsl_matrix_complex_get(A,i,j).dat[1] << " ";
        }
        std::cout << std::endl;
    }    
    */
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            std::cout << gsl_matrix_complex_get(C,i,j).dat[0] << " " << gsl_matrix_complex_get(A,i,j).dat[0] << std::endl;
            std::cout << gsl_matrix_complex_get(C,i,j).dat[1] << " " << gsl_matrix_complex_get(A,i,j).dat[1] << std::endl;
        }
    }

    return 0;
}