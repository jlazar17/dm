#include "neuosc.h"

// Macros
#define SQR(x)      ((x)*(x))                        // x^2
#define SQR_ABS(x)  (SQR(creal(x)) + SQR(cimag(x)))  // |x|^2
#define POW10(x)    (exp(M_LN10*(x)))                // 10^x
#define MIN(X,Y)    ( ((X) < (Y)) ? (X) : (Y) )
#define MAX(X,Y)    ( ((X) > (Y)) ? (X) : (Y) )
#define SIGN(a,b)   ( (b) > 0.0 ? (fabs(a)) : (-fabs(a)) )
#define KRONECKER(i,j)  ( (i)==(j) ? 1 : 0 )

gsl_matrix_complex* R(int i,int j, int cp, PhysConst* param){
    int k;
    double sd,cd,thcp;
    gsl_complex faseCP;
    gsl_complex sij,cij;
    
    gsl_matrix_complex *R = gsl_matrix_complex_calloc(param->numneu,param->numneu);
    gsl_matrix_complex_set_zero(R);
    
    // no funny business - strict order
    if (j < i) {
        k = i;
        i = j;
        j = k;
        }
        
    // reading values from params
    sij = gsl_complex_rect(gsl_matrix_get(param->s,i,j),0.0);
    cij = gsl_complex_rect(gsl_matrix_get(param->c,i,j),0.0);
    
    // diagonal terms
    for(k=0;k<param->numneu;k++){
        if (k != i-1 && k != j-1) {
            gsl_matrix_complex_set(R,k,k,gsl_complex_rect(1.0,0.0));
            }
        else {
            gsl_matrix_complex_set(R,k,k,cij);
            }
        };
    // non-diagonal terms
    if (cp != 0 ){
        thcp = gsl_matrix_get(param->dcp,cp,0);
        sd = sin(thcp);
        cd = cos(thcp);
        faseCP = gsl_complex_rect(cd,sd);
        }
    else {
        faseCP = gsl_complex_rect(1.0,0.0);
    };
    
    gsl_matrix_complex_set(R,i-1,j-1,gsl_complex_mul(gsl_complex_mul(gsl_complex_rect(1.0,0.0),sij),gsl_complex_conjugate(faseCP)));
    gsl_matrix_complex_set(R,j-1,i-1,gsl_complex_mul(gsl_complex_mul(gsl_complex_rect(-1.0,0.0),sij),faseCP));

    return R;
    }

gsl_matrix_complex* MixMatrix(PhysConst* param){
    gsl_matrix_complex *U = gsl_matrix_complex_alloc(param->numneu,param->numneu);
    if (param->numneu == 2){
        U = R(1,2,0,param);
        
        #ifdef MixMatrix_DEBUG
        printf("MixMatrix:DEBUG:R12 \n");
        gsl_matrix_complex_fprintf(stdout,R12,"%g");
        #endif
    } else if (param->numneu == 3){
        // calculating rotations
        gsl_matrix_complex* R12 = R(1,2,0,param);
        gsl_matrix_complex* R13 = R(1,3,1,param);
        gsl_matrix_complex* R23 = R(2,3,0,param);
        // tmp matrix
        gsl_matrix_complex *T1 = gsl_matrix_complex_alloc(param->numneu,param->numneu);
        
        #ifdef MixMatrix_DEBUG
        printf("MixMatrix:DEBUG:R12 \n");
        gsl_matrix_complex_fprintf(stdout,R12,"%g");
        printf("MixMatrix:DEBUG:R13 \n");
        gsl_matrix_complex_fprintf(stdout,R13,"%g");
        printf("MixMatrix:DEBUG:R23 \n");
        gsl_matrix_complex_fprintf(stdout,R23,"%g");        
        #endif
        
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R13,R12,gsl_complex_rect(0.0,0.0),T1);
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R23,T1,gsl_complex_rect(0.0,0.0),U);
        
        gsl_matrix_complex_free(R12);
        gsl_matrix_complex_free(R13);
        gsl_matrix_complex_free(R23);
        gsl_matrix_complex_free(T1);        
    } else if (param->numneu == 4) {
        // calculating rotations
        gsl_matrix_complex* R12 = R(1,2,0,param);
        gsl_matrix_complex* R13 = R(1,3,1,param);
        gsl_matrix_complex* R23 = R(2,3,0,param);
        gsl_matrix_complex* R14 = R(1,4,0,param);
        gsl_matrix_complex* R24 = R(2,4,2,param);
        gsl_matrix_complex* R34 = R(3,4,0,param);        
        // tmp matrix
        gsl_matrix_complex *T1 = gsl_matrix_complex_alloc(param->numneu,param->numneu);
        gsl_matrix_complex *T2 = gsl_matrix_complex_alloc(param->numneu,param->numneu);
        
        // multiplying matrix
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R13,R12,gsl_complex_rect(0.0,0.0),T1);
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R23,T1,gsl_complex_rect(0.0,0.0),T2);
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R14,T2,gsl_complex_rect(0.0,0.0),T1);
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R24,T1,gsl_complex_rect(0.0,0.0),T2);
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R34,T2,gsl_complex_rect(0.0,0.0),U);
        
        gsl_matrix_complex_free(R12);
        gsl_matrix_complex_free(R13);
        gsl_matrix_complex_free(R23);
        gsl_matrix_complex_free(R14);
        gsl_matrix_complex_free(R24);
        gsl_matrix_complex_free(R34);
        gsl_matrix_complex_free(T1);
        gsl_matrix_complex_free(T2);
        
    } else if (param->numneu == 5) {
        // calculating rotations
        gsl_matrix_complex* R12 = R(1,2,0,param);
        gsl_matrix_complex* R13 = R(1,3,1,param);
        gsl_matrix_complex* R23 = R(2,3,0,param);
        gsl_matrix_complex* R14 = R(1,4,2,param);
        gsl_matrix_complex* R24 = R(2,4,0,param);
        gsl_matrix_complex* R34 = R(3,4,0,param);
        gsl_matrix_complex* R15 = R(1,5,3,param);
        gsl_matrix_complex* R25 = R(2,5,0,param);
        gsl_matrix_complex* R35 = R(3,5,0,param);
        gsl_matrix_complex* R45 = R(4,5,0,param);
        // tmp matrix
        gsl_matrix_complex *T1 = gsl_matrix_complex_alloc(param->numneu,param->numneu);
        gsl_matrix_complex *T2 = gsl_matrix_complex_alloc(param->numneu,param->numneu);
               
        // multiplying matrix
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R13,R12,gsl_complex_rect(0.0,0.0),T1);
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R23,T1,gsl_complex_rect(0.0,0.0),T2);
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R14,T2,gsl_complex_rect(0.0,0.0),T1);
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R24,T1,gsl_complex_rect(0.0,0.0),T2);
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R34,T2,gsl_complex_rect(0.0,0.0),T1);
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R15,T1,gsl_complex_rect(0.0,0.0),T2);
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R25,T2,gsl_complex_rect(0.0,0.0),T1);        
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R35,T1,gsl_complex_rect(0.0,0.0),T2);
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R45,T2,gsl_complex_rect(0.0,0.0),U);
        
        gsl_matrix_complex_free(R12);
        gsl_matrix_complex_free(R13);
        gsl_matrix_complex_free(R23);
        gsl_matrix_complex_free(R14);
        gsl_matrix_complex_free(R24);
        gsl_matrix_complex_free(R34);
        gsl_matrix_complex_free(R15);
        gsl_matrix_complex_free(R25);
        gsl_matrix_complex_free(R35);
        gsl_matrix_complex_free(R45);
        gsl_matrix_complex_free(T1);
        gsl_matrix_complex_free(T2);             

    } else if (param->numneu == 6) {
        // calculating rotations
        gsl_matrix_complex* R12 = R(1,2,0,param);
        gsl_matrix_complex* R13 = R(1,3,1,param);
        gsl_matrix_complex* R23 = R(2,3,0,param);
        gsl_matrix_complex* R14 = R(1,4,0,param);
        gsl_matrix_complex* R25 = R(2,5,0,param);
        gsl_matrix_complex* R36 = R(3,6,0,param);        
        // tmp matrix
        gsl_matrix_complex *T1 = gsl_matrix_complex_alloc(param->numneu,param->numneu);
        gsl_matrix_complex *T2 = gsl_matrix_complex_alloc(param->numneu,param->numneu);                
        // multiplying matrix
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R13,R12,gsl_complex_rect(0.0,0.0),T1);
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R23,T1,gsl_complex_rect(0.0,0.0),T2);
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R14,T2,gsl_complex_rect(0.0,0.0),T1);
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R25,T1,gsl_complex_rect(0.0,0.0),T2);
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R36,T2,gsl_complex_rect(0.0,0.0),U);
        
        gsl_matrix_complex_free(R12);
        gsl_matrix_complex_free(R13);
        gsl_matrix_complex_free(R23);
        gsl_matrix_complex_free(R14);
        gsl_matrix_complex_free(R25);
        gsl_matrix_complex_free(R36);
        gsl_matrix_complex_free(T1);
        gsl_matrix_complex_free(T2);
    } else {
        cout << "Sorry, too many neutrinos. Not yet implemented! =(." << endl;
        exit(0);
    }

    if (param->neutype == 1){
        gsl_matrix_complex_conjugate(U);
    }
        
    return U;
}

gsl_matrix_complex* massM2(PhysConst* param){
    gsl_matrix_complex *M2 = gsl_matrix_complex_alloc(param->numneu,param->numneu);
    gsl_matrix_complex_set_zero(M2);    
    gsl_complex dmsq;

    for(int k=1;k<param->numneu;k++){
        dmsq = gsl_complex_rect(gsl_matrix_get(param->dmsq,k,0),0.0);
        gsl_matrix_complex_set(M2,k,k,dmsq);
    }
    return M2;
}

gsl_matrix_complex* flavorM2(PhysConst* param){
    int numneu = param->numneu;
    
    gsl_matrix_complex *U = MixMatrix(param);
    gsl_matrix_complex *U_copy = gsl_matrix_complex_alloc(numneu,numneu);
    gsl_matrix_complex_memcpy(U_copy,U);
    gsl_matrix_complex *mass2 = massM2(param);
    
    gsl_matrix_complex *T1 = gsl_matrix_complex_alloc(numneu,numneu);
    gsl_matrix_complex *T2 = gsl_matrix_complex_alloc(numneu,numneu);
    
    #ifdef flavorM2_DEBUG
    printf("flavorM2:DEBUG:Mixing Matrix \n");
    gsl_matrix_complex_fprintf(stdout,U,"%g");
    printf("flavorM2:DEBUG:Square mass diferences [mass basis] \n");
    gsl_matrix_complex_fprintf(stdout,mass2,"%g"); 
    #endif
    
    gsl_blas_zgemm(CblasNoTrans,CblasConjTrans,
                   gsl_complex_rect(1.0,0.0),mass2,U,
                   gsl_complex_rect(0.0,0.0),T1);
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,
                   gsl_complex_rect(1.0,0.0),U_copy,
                   T1,gsl_complex_rect(0.0,0.0),T2);
    
    gsl_matrix_complex_free(mass2);
    gsl_matrix_complex_free(T1);
    gsl_matrix_complex_free(U);
    gsl_matrix_complex_free(U_copy);
    
    return T2;
}

gsl_matrix_complex* flavorAcc(PhysConst* param,double E,Body *body,Track *track){
    gsl_matrix_complex *A = gsl_matrix_complex_alloc(param->numneu,param->numneu);
    gsl_matrix_complex_set_zero(A);
    
    double ye = body->ye(track);
    double density = body->density(track);
    
    #ifdef flavorAcc_DEBUG
    printf("flacorAcc:DEBUG:Body position : %g \n",track->x);
    #endif
    
    #ifdef flavorAcc_DEBUG
    printf("flavorAcc:DEBUG:Body data : \n");
    printf(" density : %g \n",density);
    printf(" ye : %g \n",ye);
    printf("conversion factor : %g \n", param->sqrt2*param->GF*param->Na*pow(param->cm,-3));
    #endif
    
    double CC = param->sqrt2*param->GF*param->Na*pow(param->cm,-3)*density*ye;
    double NC = CC*(-0.5*(1.0-ye)/ye);
    
    #ifdef flavorAcc_DEBUG
    printf("flavorAcc:DEBUG:Potentials : \n");
    printf(" CC: %g \n",CC);
    printf(" NC: %g \n",NC);
    #endif
    
    if(param->neutype == 0){
        gsl_matrix_complex_set(A,0,0,gsl_complex_rect(CC+NC,0.0));
        gsl_matrix_complex_set(A,1,1,gsl_complex_rect(NC,0.0));
        if(param->numneu>2){
            gsl_matrix_complex_set(A,2,2,gsl_complex_rect(NC,0.0));    
        }
    } else {
        gsl_matrix_complex_set(A,0,0,gsl_complex_rect(-CC-NC,0.0));
        gsl_matrix_complex_set(A,1,1,gsl_complex_rect(-NC,0.0));
        if(param->numneu>2){
            gsl_matrix_complex_set(A,2,2,gsl_complex_rect(-NC,0.0));    
        }
    }
    return A;
}

gsl_matrix_complex* massAcc(PhysConst* param,double E,Body *body,Track *track,gsl_matrix_complex* U){
    int numneu = param->numneu;
    gsl_matrix_complex *A = flavorAcc(param,E,body,track);
    // U is the mixing matrix that writes the flavor basis a function of the mass basis
    gsl_matrix_complex *U1 = gsl_matrix_complex_alloc(numneu,numneu);
    gsl_matrix_complex *U2 = gsl_matrix_complex_alloc(numneu,numneu);
    gsl_matrix_complex_memcpy(U1,U);
    gsl_matrix_complex_memcpy(U2,U);    
    // temporal matrices to make the multiplication
    gsl_matrix_complex *T1 = gsl_matrix_complex_alloc(numneu,numneu);
        
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,
                   gsl_complex_rect(1.0,0.0),A,
                   U1,gsl_complex_rect(0.0,0.0),T1);
    gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,
                   gsl_complex_rect(1.0,0.0),U2,
                   T1,gsl_complex_rect(0.0,0.0),A);
    
    gsl_matrix_complex_free(U1);
    gsl_matrix_complex_free(U2);  
    gsl_matrix_complex_free(T1);
    
    return A;
}

gsl_matrix_complex* flavorH(PhysConst* param,double E,Body *body,Track *track,gsl_matrix_complex* fM2){
    gsl_matrix_complex *I_matrix = gsl_matrix_complex_alloc(param->numneu,param->numneu);
    gsl_matrix_complex_set_identity(I_matrix);
    
    #ifdef NoOscMatterEffect
        gsl_matrix_complex *A = gsl_matrix_complex_calloc(param->numneu,param->numneu);
    #else
        gsl_matrix_complex *A = flavorAcc(param,E,body,track);
    #endif
    
    #ifdef flavorH_DEBUG
        printf("flavorH:DEBUG:Matter potential \n");
        gsl_matrix_complex_fprintf(stdout,A,"%g");
        printf("flavorH:DEBUG:Square mass diferences \n");
        gsl_matrix_complex_fprintf(stdout,fM2,"%g"); 
    #endif
        
    // multiplying matrix. result returned in A.
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0/(2.0*E),0.0),fM2,I_matrix,gsl_complex_rect(1.0,0.0),A);

    
    #ifdef Include_NewPhysics    
        double delta_NP = 1.0e-27;
        double NP_power = 1.0;
        double NP_coeff = delta_NP*pow(E,NP_power)/2.0;
        
        gsl_matrix_complex* U_NP_1 = R(param->muon+1,param->tau+1,0,param);
        gsl_matrix_complex* U_NP_2 = R(param->muon+1,param->tau+1,0,param);
            
        gsl_matrix_complex* NP_mass = gsl_matrix_complex_calloc(param->numneu,param->numneu);
            
        gsl_matrix_complex_set(NP_mass,param->muon,param->muon,gsl_complex_rect(-NP_coeff,0.0));
        gsl_matrix_complex_set(NP_mass,param->tau,param->tau,gsl_complex_rect(NP_coeff,0.0));
        
        gsl_matrix_complex *T1 = gsl_matrix_complex_calloc(param->numneu,param->numneu);
        gsl_matrix_complex *NP_flavor = gsl_matrix_complex_calloc(param->numneu,param->numneu);
        
        // calculating U*NPhysics*U_dagger
        gsl_blas_zgemm(CblasNoTrans,CblasConjTrans,gsl_complex_rect(1.0,0.0),NP_mass,U_NP_1,gsl_complex_rect(0.0,0.0),T1);
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),U_NP_2,T1,gsl_complex_rect(0.0,0.0),NP_flavor);
        
        #ifdef NewPhysics_DEBUG
            cout << "== NewPhysics Term ==" << endl;
            for(int i = 0; i < param->numneu; i++){
                for(int j = 0; j < param->numneu; j++){
                    cout << gsl_matrix_complex_get(NP_flavor,i,j).dat[0] <<
                        "+i" << gsl_matrix_complex_get(NP_flavor,i,j).dat[1] << " ";
                } cout << endl;
            } 
            cout << "==" << endl;
        #endif
        
        gsl_matrix_complex_add(A,NP_flavor);
        
        gsl_matrix_complex_free(U_NP_1);
        gsl_matrix_complex_free(U_NP_2);
        gsl_matrix_complex_free(NP_mass);
        gsl_matrix_complex_free(NP_flavor);
        gsl_matrix_complex_free(T1);
    #endif
    
    #ifdef flavorH_DEBUG
    printf("flavorH:DEBUG:Hamiltonian \n");
    gsl_matrix_complex_fprintf(stdout,A,"%g");
    #endif    
    
    gsl_matrix_complex_free(I_matrix);
    
    return A;
}

gsl_matrix_complex* massH(PhysConst* param,double E,Body *body,Track *track,gsl_matrix_complex* mM2,gsl_matrix_complex* U){
    gsl_matrix_complex *I_matrix = gsl_matrix_complex_alloc(param->numneu,param->numneu);
    gsl_matrix_complex_set_identity(I_matrix);
    
    #ifdef NoOscMatterEffect
        gsl_matrix_complex *A = gsl_matrix_complex_calloc(param->numneu,param->numneu);
    #else
        gsl_matrix_complex *A = massAcc(param,E,body,track,U);
    #endif
    
    #ifdef flavorH_DEBUG
        printf("flavorH:DEBUG:Matter potential \n");
        gsl_matrix_complex_fprintf(stdout,A,"%g");
        printf("flavorH:DEBUG:Square mass diferences \n");
        gsl_matrix_complex_fprintf(stdout,fM2,"%g"); 
    #endif
        
    // multiplying matrix. result returned in A.
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0/(2.0*E),0.0),mM2,I_matrix,gsl_complex_rect(1.0,0.0),A);
    
    gsl_matrix_complex_free(I_matrix);
    
    return A;
}

gsl_matrix_complex* massExpH0(PhysConst* param,double E,Track *track,gsl_matrix_complex* massM2){
    // this rountine calculates the exponential
    // of the constant part of the Hamiltonian
    // in the mass basis, namely,
    // EXP(i (\Delta m^2 / 2E) * x ) 
    gsl_matrix_complex *exp_matrix = gsl_matrix_complex_calloc(param->numneu,param->numneu);
    
    gsl_complex exp_term;

    for(int k=0;k<param->numneu;k++){
        exp_term = gsl_complex_polar(1.0,
                               exp(gsl_matrix_get(param->dmsq,k,0)*track->x/(2.0*E)));
        gsl_matrix_complex_set(exp_matrix,k,k,exp_term);
    }
    
    //gsl_matrix_complex_print(exp_matrix);
    
    return exp_matrix;
}

/*
int RHS_INT_GSL(double x,gsl_complex_vector neuvec,ParamsContainer params){
    double x_prev = params.rpos;
    double x_current = x;
    params.track.x = x_current;
    
    //gsl_matrix_complex_memcpy
    
    gsl_matrix_complex *H_current = flavorH(params.param,params.E,params.body,params.track,params.flavorM2);
    gsl_matrix_complex *H1 = H_current = params.H0;
    
    // calculating the eigen values and eigen vectors of H1
    // save copying H1
    gsl_matrix_complex H1
    gsl_matrix_complex_memcpy(H_copy,H1);
    // allocating space for the calculation
    gsl_eigen_herm_workspace *w = gsl_eigen_herm_alloc(param.numneu);
    // calculating hermitian matrix
    gsl_eigen_hermv(H1,lambda,Q,w);
    // freeing allocated memspace
    gsl_eigen_herm_free(w);
    
    gsl_matrix_complex *S = ;
    gsl_matrix_complex * interaction_H = ; 
    
    rhs_int = ;
    
    rhs = rhs_int,neuvec;
    
    params.cH = H_current;
    params.cS = S;
    
    return rhs;
}
*/

int RHS_GSL(double x,const double *neuvec_real_in, double *neuvec_real_out,void *params_void){
    // typing params
    Container *params = (Container *) params_void;
    // dimentions
    int numneu = params->param->numneu;
    // initializing complex vectors
    gsl_vector_complex *neuvec_complex_in = gsl_vector_complex_calloc (numneu);
    gsl_vector_complex *neuvec_complex_out = gsl_vector_complex_calloc (numneu);
    // copying 'real vector' to 'complex vector'
    for(int k=0;k<numneu;k++){
        gsl_vector_complex_set(neuvec_complex_in,k,gsl_complex_rect(neuvec_real_in[2*k],neuvec_real_in[2*k+1]));
    }
    
    #ifdef RHS_GSL_DEBUG
    printf("RHS_GSL:DEBUG:Real entry neutrino state (");
    for (int j = 0; j < numneu; j++)
        printf(" %g+%g*i ", neuvec_real_in[2*j],neuvec_real_in[2*j+1]);
    printf(")\n");
    printf("RHS_GSL:DEBUG:Complex entry neutrino state \n");
    gsl_vector_complex_fprintf(stdout,neuvec_complex_in,"%g");
    #endif
    
    // updating position
    params->track->x = x;
    // calculating hamiltoninan
    gsl_matrix_complex *H_current = flavorH(params->param,params->E,params->body,params->track,params->flavorM2);
    // multiplying    
    gsl_complex ii = gsl_complex_rect(0.0,-1.0);
    gsl_complex zero = gsl_complex_rect(0.0,0.0);
    gsl_blas_zgemv(CblasNoTrans,ii,H_current,neuvec_complex_in,zero,neuvec_complex_out);
    
    #ifdef RHS_GSL_DEBUG
    printf("RHS_GSL:DEBUG:Current Hamiltonian \n");
    gsl_matrix_complex_fprintf(stdout,H_current,"%g");
    #endif
    
    // copying 'complex vector' to 'real vector'
    gsl_complex component;
    for(int k=0;k<numneu;k++){
        component = gsl_vector_complex_get (neuvec_complex_out,k);
        neuvec_real_out[2*k] = GSL_REAL(component);
        neuvec_real_out[2*k+1] = GSL_IMAG(component);
    }
    
    #ifdef RHS_GSL_DEBUG
    printf("RHS_GSL:DEBUG:Complex exit neutrino state (");
    for (int j = 0; j < numneu; j++)
        printf(" %g+%g*i ", neuvec_real_out[2*j],neuvec_real_out[2*j+1]);
    printf(")\n");
    printf("RHS_GSL:DEBUG:Complex exit neutrino state \n");
    gsl_vector_complex_fprintf(stdout,neuvec_complex_out,"%g");
    #endif
    
    gsl_matrix_complex_free(H_current);
    gsl_vector_complex_free(neuvec_complex_in);
    gsl_vector_complex_free(neuvec_complex_out);
    
    return GSL_SUCCESS;
}

int CalNeuOscGSL(double osc_prob[],int ineu,double E,Track *track,Body *body,gsl_matrix_complex* flavorM2,PhysConst* param,double abs_error,double rel_error,bool return_state = false,bool optimization = true){
    int gsl_status = GSL_SUCCESS;
    int numneu = param->numneu;
    
    // creating auxiliary variable    
    Container* params = new Container;
    params->track = track;
    params->body = body;
    params->param = param;
    params->flavorM2 = flavorM2;
    //params.flavorM2 = gsl_matrix_complex_alloc(numneu,numneu);
    //gsl_matrix_complex_memcpy(params.flavorM2,flavorM2);
    params->E = E;

#ifdef CalNeuOscGSL_DEBUG
    printf("Oscillation paramers :\n");
    printf("E : %g [eV]\n", params->E);
    printf("xini : %g [eV^-1]\n", params->track->xini);
    printf("xend : %g [eV^-1]\n", params->track->xend);
    printf("Square mass difference matrix \n");
    gsl_matrix_complex_fprintf(stdout,params->flavorM2,"%g");
    
    printf("Get body information \n");
    cout << "name : " << params->body->name << endl;
    printf("density : %g [gr/cm^3]\n",params->body->density(params->track));
#endif
    
    // setting up GSL ODE solver
    gsl_odeiv_step *s = gsl_odeiv_step_alloc(gsl_odeiv_step_rkf45, 2*numneu);
    gsl_odeiv_control *c = gsl_odeiv_control_y_new(rel_error,abs_error);
    gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc(2*numneu);
    gsl_odeiv_system sys = {&RHS_GSL, NULL, 2*numneu,params};
    
    //gsl_odeiv_step_reset(s);
    //gsl_odeiv_evolve_reset(e);
    
    // initial neutrino state in flavor space
    double neuvec_real[2*numneu];
    for (int j=0; j< 2*numneu; j++)
        neuvec_real[j] = 0.0;
    neuvec_real[2*ineu] = 1.0;

#ifdef CalNeuOscGSL_DEBUG
    printf("Initial neutrino state: (");
    for (int j = 0; j < numneu; j++)
        printf(" %g+%g*i ", neuvec_real[2*j],neuvec_real[2*j+1]);
    printf(")\n");
#endif
    
    // defining ODE extra variables
    double x = 0;               // ODE independent variable
    double x_ini = params->track->xini;  // initial position
    double x_end = params->track->xend;  // final position
    // step sizes
    double h        = MIN(1.0*param->km,x_end/10.0);
    double h_min    = 1.0e-5*param->km;
    
#ifdef CalNeuOscGSL_DEBUG
    printf("GSL paramers :\n");
    printf("x_ini : %g [eV^-1]\n", x_ini);
    printf("x_end : %g [eV^-1]\n", x_end);
    printf("h : %g [eV^-1]\n", h);
    printf("h_min : %g [eV^-1]\n", h_min);        
#endif
    // initial position
    x = x_ini;
    
    while (x < x_end){

    gsl_status = gsl_odeiv_evolve_apply(e,c,s,&sys,&x,x_end,&h,neuvec_real);

#ifdef CalNeuOscGSLStep_DEBUG
    printf("x_current : %g [eV^-1]\n", x);
    printf("step : %g [eV^-1]\n", h);    
#endif 
     
    if( gsl_status != GSL_SUCCESS ){
        fprintf(stderr,"CalNeuOscGSL: Error in GSL ODE solver,\n");
        break;
    }
    if ( h - h_min < 0.0 ){
        h = h_min;
    }
    }
    
#ifdef CalNeuOscGSL_DEBUG
    printf("Final neutrino state: (");
    for (int j = 0; j < numneu; j++)
        printf(" %g+%g*i ", neuvec_real[2*j],neuvec_real[2*j+1]);
    printf(")\n");
#endif    
    
    if (e) { gsl_odeiv_evolve_free(e);      e = NULL; }
    if (c) { gsl_odeiv_control_free(c);     c = NULL; }
    if (s) { gsl_odeiv_step_free(s);        s = NULL; }
    
    if (return_state) {
        for (int j=0; j< 2*numneu; j++){
            osc_prob[j] = neuvec_real[j];
        }
    } else {
        for (int j=0; j< numneu; j++){
            osc_prob[j] = SQR(neuvec_real[2*j])+SQR(neuvec_real[2*j+1]);
        }
    }
    
    free(params);
    
    return GSL_SUCCESS;
}
