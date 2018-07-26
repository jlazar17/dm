#include "gpu_tools.h"
#define SQR(x)      ((x)*(x))                        // x^2

int CUDAComplexMatrixAdd(gsl_matrix_complex* A_gsl,gsl_matrix_complex* B_gsl){
    // cublas variables
    cudaError_t cuda_status;
    cublasStatus_t status;
    cublasHandle_t handle;

    //cuDoubleComplex* A_cuda;
    double* A_cuda;
    double* B_cuda;
    
    int complex_matrix_size = A_gsl->size1;
    int vector_size = 2*SQR(complex_matrix_size);
    
    // get complex data
    double* A_dbl = A_gsl->data;
    double* B_dbl = B_gsl->data;
    
    double alpha = 1.0;
    
    /* Allocate device memory for matrices */
    cuda_status = cudaMalloc ((void**)&A_cuda, vector_size*sizeof(double));
    if (cuda_status != cudaSuccess) {
        printf ("device memory allocation failed \n");
        return EXIT_FAILURE;
    }

    cuda_status = cudaMalloc ((void**)&B_cuda, vector_size*sizeof(double));
    if (cuda_status != cudaSuccess) {
        printf ("device memory allocation failed \n");
        return EXIT_FAILURE;
    }
 
    /* Initialize CUBLAS */
    status = cublasCreate(&handle);
    if (status != CUBLAS_STATUS_SUCCESS) {
        printf ("CUBLAS initialization failed\n");
        return EXIT_FAILURE;
    }

    /* Copy matrix from host to device */
    
    status = cublasSetVector (vector_size, sizeof(double), A_dbl, 1, A_cuda, 1);
    if (status != CUBLAS_STATUS_SUCCESS) {
        printf ("data download failed\n");
        cublasDestroy(handle);
        return EXIT_FAILURE;
    }
    
    status = cublasSetVector (vector_size, sizeof(double), B_dbl, 1, B_cuda, 1);
    if (status != CUBLAS_STATUS_SUCCESS) {
        printf ("data download failed\n");
        cublasDestroy(handle);
        return EXIT_FAILURE;
    }

    /* Performs operation using cublas */
    
    status = cublasDaxpy(handle,
                          vector_size,
                          &alpha,
                          A_cuda, 1,
                          B_cuda, 1);
    
    if (status != CUBLAS_STATUS_SUCCESS) {
        printf("cublasZgeam execution error.\n");
        return EXIT_FAILURE;
    }
            
    /* Reading back the result from device memory */
    status = cublasGetVector (vector_size, sizeof(double), B_cuda, 1, A_dbl, 1);
    if (status != CUBLAS_STATUS_SUCCESS) {
        printf ("data upload failed");
        cublasDestroy(handle);        
        return EXIT_FAILURE;
    }    
    
    /* free memory */
    cudaFree(A_cuda);
    cudaFree(B_cuda);
    cublasDestroy(handle);
    
    #ifdef CUDAComplexMatrixAdd_DEBUG
    for(int i = 0; i < A_gsl->size1; i++){
        for(int j = 0; j < A_gsl->size1; j++){
            cout << gsl_matrix_complex_get(A_gsl,i,j).dat[0] <<
            "+i" << gsl_matrix_complex_get(A_gsl,i,j).dat[1] << " ";
        }
        cout << endl;
    }    
    #endif
    
    return EXIT_SUCCESS;
};


// this first attempt used matrix operation

int CUDAComplexMatrixAdd_old(gsl_matrix_complex* A_gsl,gsl_matrix_complex* B_gsl){
    // cublas variables
    cudaError_t cuda_status;
    cublasStatus_t status;
    cublasHandle_t handle;

    //cuDoubleComplex* A_cuda;
    double* A_cuda;
    double* B_cuda;
    double* C_cuda;
    // assuming square matrices converting double (2) to double
    //cout << A_gsl->size1 << " " << A_gsl->size2 << " " << A_gsl->tda << endl;
    
    int complex_matrix_size = A_gsl->size1;
    int col_size = complex_matrix_size;
    int row_size = complex_matrix_size;
    
    // get complex data
    double* A_dbl = A_gsl->data;
    double* B_dbl = B_gsl->data;
    //double* C_dbl = (double*) malloc(SQR(matrix_size)*sizeof(double));
    
    //cuDoubleComplex *alpha = make_cuDoubleComplex(1.0,0.0);
    double alpha = 1.0;
    
    /* Allocate device memory for matrices */
    cuda_status = cudaMalloc ((void**)&A_cuda, col_size*row_size*sizeof(double));
    if (cuda_status != cudaSuccess) {
        printf ("device memory allocation failed \n");
        return EXIT_FAILURE;
    }

    cuda_status = cudaMalloc ((void**)&B_cuda, col_size*row_size*sizeof(double));
    if (cuda_status != cudaSuccess) {
        printf ("device memory allocation failed \n");
        return EXIT_FAILURE;
    }

    cuda_status = cudaMalloc ((void**)&C_cuda, col_size*row_size*sizeof(double));
    if (cuda_status != cudaSuccess) {
        printf ("device memory allocation failed \n");
        return EXIT_FAILURE;
    }
    
    /* Initialize CUBLAS */
    status = cublasCreate(&handle);
    if (status != CUBLAS_STATUS_SUCCESS) {
        printf ("CUBLAS initialization failed\n");
        return EXIT_FAILURE;
    }

    /* Copy matrix from host to device */
    
    status = cublasSetMatrix (row_size, col_size, sizeof(double), A_dbl, col_size, A_cuda, col_size);
    if (status != CUBLAS_STATUS_SUCCESS) {
        printf ("data download failed\n");
        cublasDestroy(handle);
        return EXIT_FAILURE;
    }
    
    status = cublasSetMatrix (row_size, col_size, sizeof(double), B_dbl, col_size, B_cuda, col_size);
    if (status != CUBLAS_STATUS_SUCCESS) {
        printf ("data download failed\n");
        cublasDestroy(handle);
        return EXIT_FAILURE;
    }

    /* Performs operation using cublas */
    
    status = cublasDgeam(handle,
                          CUBLAS_OP_N, CUBLAS_OP_N,
                          row_size, col_size,
                          &alpha,
                          A_cuda, col_size,
                          &alpha,
                          B_cuda, col_size,
                          C_cuda, col_size);
    
    if (status != CUBLAS_STATUS_SUCCESS) {
        printf("cublasZgeam execution error.\n");
        return EXIT_FAILURE;
    }
            
    /* Reading back the result from device memory */
    status = cublasGetMatrix (row_size, col_size, sizeof(double), C_cuda, col_size, A_dbl, col_size);
    if (status != CUBLAS_STATUS_SUCCESS) {
        printf ("data upload failed");
        cublasDestroy(handle);        
        return EXIT_FAILURE;
    }    
    
    /* free memory */
    cudaFree(A_cuda);
    cudaFree(B_cuda);
    cudaFree(C_cuda);
    cublasDestroy(handle);
    
    #ifdef CUDAComplexMatrixAdd_DEBUG
    for(int i = 0; i < A_gsl->size1; i++){
        for(int j = 0; j < A_gsl->size1; j++){
            cout << gsl_matrix_complex_get(A_gsl,i,j).dat[0] <<
            "+i" << gsl_matrix_complex_get(A_gsl,i,j).dat[1] << " ";
        }
        cout << endl;
    }    
    #endif
    
    return EXIT_SUCCESS;
};

