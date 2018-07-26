cimport cython
from cython_gsl cimport *
from libcpp cimport bool
from libc.stdlib cimport malloc

# my classes
cimport gslphysconst
from gslphysconst cimport PhysConst
cimport gslbody
from gslbody cimport Body,Track

# lets define an auxiliary container

cdef class Container:
    # basic storage
    cdef Track track
    cdef PhysConst param
    cdef gsl_matrix_complex *flavorM2
    cdef Body body
    cdef double E
    
    # used for optimization
    cdef gsl_matrix_complex H0
    cdef gsl_matrix_complex S
    cdef double rpos
    cdef gsl_matrix_complex cH
    cdef gsl_matrix_complex cS

# mixing matrix implementation

cdef gsl_matrix_complex* R(int i,int j, int cp, PhysConst param):
    cdef int k
    cdef double sd,cd,thcp
    cdef gsl_complex faseCP
    cdef gsl_complex sij,cij
    
    cdef gsl_matrix_complex *R = gsl_matrix_complex_calloc(param.numneu,param.numneu)
    gsl_matrix_complex_set_zero(R)
    
    ## no funny business - strict order
    if j < i :
        k = i
        i = j
        j = k
        
    ## reading values from params
    sij = gsl_complex_rect(gsl_matrix_get(param.s,i,j),0.0)
    cij = gsl_complex_rect(gsl_matrix_get(param.c,i,j),0.0)
    
    ## diagonal terms
    for k in range(param.numneu):
        if  k != i-1 and k != j-1:
            gsl_matrix_complex_set(R,k,k,gsl_complex_rect(1.0,0.0))
        else :
            gsl_matrix_complex_set(R,k,k,cij)

    ## non-diagonal terms
    if cp != 0 :
        thcp = gsl_matrix_get(param.dcp,cp,0)
        sd = sin(thcp)
        cd = cos(thcp)
        faseCP = gsl_complex_rect(cd,sd)
    else :
        faseCP = gsl_complex_rect(1.0,0.0)
    
    gsl_matrix_complex_set(R,i-1,j-1,gsl_complex_mul(gsl_complex_mul(gsl_complex_rect(1.0,0.0),sij),gsl_complex_conjugate(faseCP)))
    gsl_matrix_complex_set(R,j-1,i-1,gsl_complex_mul(gsl_complex_mul(gsl_complex_rect(-1.0,0.0),sij),faseCP))

    return R

cdef gsl_matrix_complex* MixMatrix(PhysConst param):
    cdef gsl_matrix_complex *U = gsl_matrix_complex_alloc(param.numneu,param.numneu)
    
    # define rotation matrices
    cdef gsl_matrix_complex *R12
    cdef gsl_matrix_complex *R13
    cdef gsl_matrix_complex *R23
    cdef gsl_matrix_complex *R14
    cdef gsl_matrix_complex *R24
    cdef gsl_matrix_complex *R34
    cdef gsl_matrix_complex *R15
    cdef gsl_matrix_complex *R25
    cdef gsl_matrix_complex *R35
    cdef gsl_matrix_complex *R45
    cdef gsl_matrix_complex *R36
    
    # define auxiliary matrices
    cdef gsl_matrix_complex *T1
    cdef gsl_matrix_complex *T2
    cdef gsl_matrix_complex *T3
    cdef gsl_matrix_complex *T4
    cdef gsl_matrix_complex *T5
    cdef gsl_matrix_complex *T6
    cdef gsl_matrix_complex *T7
    cdef gsl_matrix_complex *T8
    
    if param.numneu == 3:
        # calculating rotations
        R12 = R(1,2,0,param)
        R13 = R(1,3,1,param)
        R23 = R(2,3,0,param)
        # tmp matrix
        T1 = gsl_matrix_complex_alloc(param.numneu,param.numneu)
        
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R13,R12,gsl_complex_rect(0.0,0.0),T1)
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R23,T1,gsl_complex_rect(0.0,0.0),U)
        
        # free space
        gsl_matrix_complex_free(R12)
        gsl_matrix_complex_free(R13)
        gsl_matrix_complex_free(R23)
        gsl_matrix_complex_free(T1)
        
        return U   
    elif param.numneu == 4:
        # calculating rotations
        R12 = R(1,2,0,param)
        R13 = R(1,3,1,param)
        R23 = R(2,3,0,param)
        R14 = R(2,3,0,param)
        R24 = R(2,4,2,param)
        R34 = R(3,4,0,param)        
        # tmp matrix
        T1 = gsl_matrix_complex_alloc(param.numneu,param.numneu)
        T2 = gsl_matrix_complex_alloc(param.numneu,param.numneu)
        T3 = gsl_matrix_complex_alloc(param.numneu,param.numneu)
        T4 = gsl_matrix_complex_alloc(param.numneu,param.numneu)
        # multiplying matrix
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R13,R12,gsl_complex_rect(0.0,0.0),T1)
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R23,T1,gsl_complex_rect(0.0,0.0),T2)
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R14,T2,gsl_complex_rect(0.0,0.0),T3)
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R24,T3,gsl_complex_rect(0.0,0.0),T4)
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R34,T4,gsl_complex_rect(0.0,0.0),U)
        
        # free space
        gsl_matrix_complex_free(R12)
        gsl_matrix_complex_free(R13)
        gsl_matrix_complex_free(R23)
        gsl_matrix_complex_free(R14)
        gsl_matrix_complex_free(R24)
        gsl_matrix_complex_free(R34)
        gsl_matrix_complex_free(T1)
        gsl_matrix_complex_free(T2)
        gsl_matrix_complex_free(T3)
        gsl_matrix_complex_free(T4)
                
        return U
    elif param.numneu == 5 :
        # calculating rotations
        R12 = R(1,2,0,param)
        R13 = R(1,3,1,param)
        R23 = R(2,3,0,param)
        R14 = R(1,4,2,param)
        R24 = R(2,4,0,param)
        R34 = R(3,4,0,param)
        R15 = R(1,5,3,param)
        R25 = R(2,5,0,param)
        R35 = R(3,5,0,param)
        R45 = R(4,5,0,param)
        # tmp matrix
        T1 = gsl_matrix_complex_alloc(param.numneu,param.numneu)
        T2 = gsl_matrix_complex_alloc(param.numneu,param.numneu)
        T3 = gsl_matrix_complex_alloc(param.numneu,param.numneu)
        T4 = gsl_matrix_complex_alloc(param.numneu,param.numneu)
        T5 = gsl_matrix_complex_alloc(param.numneu,param.numneu)
        T6 = gsl_matrix_complex_alloc(param.numneu,param.numneu)
        T7 = gsl_matrix_complex_alloc(param.numneu,param.numneu)
        T8 = gsl_matrix_complex_alloc(param.numneu,param.numneu)               
        # multiplying matrix
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R13,R12,gsl_complex_rect(0.0,0.0),T1)
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R23,T1,gsl_complex_rect(0.0,0.0),T2)
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R14,T2,gsl_complex_rect(0.0,0.0),T3)
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R24,T3,gsl_complex_rect(0.0,0.0),T4)
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R34,T4,gsl_complex_rect(0.0,0.0),T5)
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R15,T5,gsl_complex_rect(0.0,0.0),T6)
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R25,T6,gsl_complex_rect(0.0,0.0),T7)        
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R35,T7,gsl_complex_rect(0.0,0.0),T8)
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R45,T8,gsl_complex_rect(0.0,0.0),U)
        # free space
        gsl_matrix_complex_free(R12)
        gsl_matrix_complex_free(R13)
        gsl_matrix_complex_free(R23)
        gsl_matrix_complex_free(R14)
        gsl_matrix_complex_free(R24)
        gsl_matrix_complex_free(R34)
        gsl_matrix_complex_free(R15)
        gsl_matrix_complex_free(R25)
        gsl_matrix_complex_free(R35)
        gsl_matrix_complex_free(R45)
        gsl_matrix_complex_free(T1)
        gsl_matrix_complex_free(T2)
        gsl_matrix_complex_free(T3)
        gsl_matrix_complex_free(T4)
        gsl_matrix_complex_free(T5)
        gsl_matrix_complex_free(T6)
        gsl_matrix_complex_free(T7)
        gsl_matrix_complex_free(T8)                
        
        return U
    elif param.numneu == 6:
        # calculating rotations
        R12 = R(1,2,0,param)
        R13 = R(1,3,1,param)
        R23 = R(2,3,0,param)
        R14 = R(2,3,0,param)
        R25 = R(2,5,0,param)
        R36 = R(3,6,0,param)        
        # tmp matrix
        T1 = gsl_matrix_complex_alloc(param.numneu,param.numneu)
        T2 = gsl_matrix_complex_alloc(param.numneu,param.numneu)
        T3 = gsl_matrix_complex_alloc(param.numneu,param.numneu)
        T4 = gsl_matrix_complex_alloc(param.numneu,param.numneu)                
        # multiplying matrix
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R13,R12,gsl_complex_rect(0.0,0.0),T1)
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R23,T1,gsl_complex_rect(0.0,0.0),T2)
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R14,T2,gsl_complex_rect(0.0,0.0),T3)
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R25,T3,gsl_complex_rect(0.0,0.0),T4)
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),R36,T4,gsl_complex_rect(0.0,0.0),U)
        # free space
        gsl_matrix_complex_free(R12)
        gsl_matrix_complex_free(R13)
        gsl_matrix_complex_free(R23)
        gsl_matrix_complex_free(R14)
        gsl_matrix_complex_free(R25)
        gsl_matrix_complex_free(R36)
        gsl_matrix_complex_free(T1)
        gsl_matrix_complex_free(T2)
        gsl_matrix_complex_free(T3)
        gsl_matrix_complex_free(T4)
        
        return U
    else:
        print "Sorry, too many neutrinos. Not yet implemented! =(."
        quit()

## DEFINING THE HAMILTONIAN ##

cdef gsl_matrix_complex* massM2(PhysConst param):
    cdef gsl_matrix_complex *M2 = gsl_matrix_complex_alloc(param.numneu,param.numneu)
    gsl_matrix_complex_set_zero(M2)  
    cdef gsl_complex dmsq

    for k in range(1,param.numneu):
        dmsq = gsl_complex_rect(gsl_matrix_get(param.dmsq,k,0),0.0)
        gsl_matrix_complex_set(M2,k,k,dmsq)
    
    return M2

cdef gsl_matrix_complex* flavorM2(PhysConst param):
    cdef int numneu = param.numneu
    
    cdef gsl_matrix_complex *U = MixMatrix(param)
    cdef gsl_matrix_complex *U_copy = gsl_matrix_complex_alloc(numneu,numneu)
    gsl_matrix_complex_memcpy(U_copy,U)
    cdef gsl_matrix_complex *mass2 = massM2(param)
    
    cdef gsl_matrix_complex *T1 = gsl_matrix_complex_alloc(numneu,numneu)
    cdef gsl_matrix_complex *T2 = gsl_matrix_complex_alloc(numneu,numneu)
    
    gsl_blas_zgemm(CblasNoTrans,CblasConjTrans,gsl_complex_rect(1.0,0.0),mass2,U,gsl_complex_rect(0.0,0.0),T1)
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),U_copy,T1,gsl_complex_rect(0.0,0.0),T2)
    
    gsl_matrix_complex_free(T1)
    
    return T2

cdef gsl_matrix_complex* flavorAcc(PhysConst param,double E,Body body,Track track):
    cdef gsl_matrix_complex *A = gsl_matrix_complex_alloc(param.numneu,param.numneu)
    gsl_matrix_complex_set_zero(A)
    
    cdef double ye = body.ye(track)
    cdef double density = body.density(track)
    
    cdef double CC = param.sqrt2*param.GF*param.Na*pow(param.cm,3)*density*ye
    cdef double NC = CC*(-0.5*(1.0-ye)/ye)
    
    if str(param.neutype) == "neutrino" : 
        gsl_matrix_complex_set(A,0,0,gsl_complex_rect(CC+NC,0.0))
        gsl_matrix_complex_set(A,1,1,gsl_complex_rect(NC,0.0))
        gsl_matrix_complex_set(A,2,2,gsl_complex_rect(NC,0.0))
    else : 
        gsl_matrix_complex_set(A,0,0,gsl_complex_rect(-CC-NC,0.0))
        gsl_matrix_complex_set(A,1,1,gsl_complex_rect(-NC,0.0))
        gsl_matrix_complex_set(A,2,2,gsl_complex_rect(-NC,0.0))
    
    return A

cdef gsl_matrix_complex* flavorH(PhysConst param,double E,Body body,Track track,gsl_matrix_complex* fM2):
    cdef gsl_matrix_complex *I_matrix = gsl_matrix_complex_alloc(param.numneu,param.numneu)
    gsl_matrix_complex_set_identity(I_matrix)
    cdef gsl_matrix_complex *A = flavorAcc(param,E,body,track)
    
    # multiplying matrix. result returned in A.
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0/(2.0*E),0.0),fM2,I_matrix,gsl_complex_rect(1.0,0.0),A)
    # free
    gsl_matrix_complex_free(I_matrix)
    
    return A

## RUNGE-KUTTA ##

cdef int RHS_GSL(double x, double neuvec_real_in[], double neuvec_real_out[],void *params_void):
    # typing params
    cdef Container params = <Container> params_void
    #cdef Container *params = <Container *> params_void
    # dimentions
    cdef int numneu = params.param.numneu
    # initializing complex vectors
    cdef gsl_vector_complex *neuvec_complex_in = gsl_vector_complex_calloc (numneu)
    cdef gsl_vector_complex *neuvec_complex_out = gsl_vector_complex_calloc (numneu)
    # copying 'real vector' to 'complex vector'
    for k in range(numneu):
        gsl_vector_complex_set(neuvec_complex_in,k,gsl_complex_rect(neuvec_real_in[2*k],neuvec_real_in[2*k+1]))
    
    # updating position
    params.track.x = x
    # calculating hamiltoninan
    cdef gsl_matrix_complex *H_current = flavorH(params.param,params.E,params.body,params.track,params.flavorM2)
    # multiplying    
    cdef gsl_complex ii = gsl_complex_rect(0.0,-1.0)
    cdef gsl_complex zero = gsl_complex_rect(0.0,0.0)
    gsl_blas_zgemv(CblasNoTrans,ii,H_current,neuvec_complex_in,zero,neuvec_complex_out)
    
    # copying 'complex vector' to 'real vector'
    cdef gsl_complex component
    for k in range(numneu):
        component = gsl_vector_complex_get (neuvec_complex_out,k)
        neuvec_real_out[2*k] = GSL_REAL(component)
        neuvec_real_out[2*k+1] = GSL_IMAG(component)
    
    return GSL_SUCCESS

cdef int CalNeuOscGSL(double osc_prob[],int ineu,double E,Track track,Body body,gsl_matrix_complex* flavorM2,PhysConst param,double abs_error,double rel_error,bool return_state = False,bool optimization = True):
    cdef int gsl_status = GSL_SUCCESS
    cdef int numneu = 3
    
    # creating auxiliary variable    
    cdef Container params = new Container
    # filling it up
    params.track = track
    PhysConst params.param = param
    params.flavorM2 = flavorM2
    params.body = body
    params.E = E
    
    # setting up GSL ODE solver
    cdef gsl_odeiv_step *s = gsl_odeiv_step_alloc(gsl_odeiv_step_rkf45, 2*numneu)
    cdef gsl_odeiv_control *c = gsl_odeiv_control_y_new(rel_error,abs_error)
    cdef gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc(2*numneu)
    cdef gsl_odeiv_system sys = {&RHS_GSL, NULL, 2*numneu,&params}
    
    #gsl_odeiv_step_reset(s);
    #gsl_odeiv_evolve_reset(e);
    
    # initial neutrino state in flavor space
    cdef int numneu2 = 2*numneu
    cdef double* neuvec_real = malloc(numneu2*cython.sizeof(double))
    for j in range(2*numneu):
        neuvec_real[j] = 0.0
    neuvec_real[2*ineu] = 1.0
    
    # defining ODE extra variables
    cdef double x = 0               # ODE independent variable
    cdef double x_ini = track.xini  # initial position
    cdef double x_end = track.xend  # final position
    # step sizes
    cdef double h        = MIN(1.0*param.km,x_end/10.0)
    cdef double h_min    = 1.0e-5*param.km

    # initial position
    x = x_ini
    
    while x < x_end :
        gsl_status = gsl_odeiv_evolve_apply(e,c,s,&sys,&x,x_end,&h,neuvec_real);
         
        if gsl_status != GSL_SUCCESS :
            print "CalNeuOscGSL: Error in GSL ODE solver"
            quit()
            
        if h < h_min :
            h = h_min
    
    if e :
        gsl_odeiv_evolve_free(e)
        e = NULL
    if c : 
        gsl_odeiv_control_free(c)
        c = NULL
    if s :
        gsl_odeiv_step_free(s)
        s = NULL
    
    if return_state :
        osc_prob = neuvec_real    
    else :
        for j in range(numneu):
            osc_prob[j] = SQR(neuvec_real[2*j])+SQR(neuvec_real[2*j+1])
    
    return GSL_SUCCESS

def main():
    cdef gsl_complex c2
    GSL_SET_COMPLEX(&c2, 2.2, 1.5)
    print 'c2',GSL_REAL(c2)