#include "SUNalg.h"
/*
 * We implement the SU(N) semialgebraic solution of the
 * diffusion problem. We wil asumme that the problem
 * will be solve in the mass-interaction basis
*/

/*
-----------------------------------------------------------------------
SU_vector implementation
-----------------------------------------------------------------------         
*/

SU_vector::SU_vector(){
    dim = 0;
    memset(components, 0.0, sizeof(double)*NEUMAXSQ);
};

SU_vector::SU_vector(int d){
    dim = d;
    //components = new double[size];
    memset(components, 0.0, sizeof(double)*NEUMAXSQ);
};

SU_vector::SU_vector(int d,double* comp){
    dim = d;
    for(int i = 0; i< SQR(dim); i++)
        components[i] = comp[i];
    //components = comp;
};

SU_vector::SU_vector(std::vector<double> comp){
    dim = (int) sqrt(comp.size());
    for(int i = 0; i< SQR(dim); i++)
        components[i] = comp[i];
    //components = comp;
};

SU_vector::SU_vector(gsl_matrix_complex* m){
    dim = m->size1;
    int size = SQR(dim);
    // create and set components to zero
    //components = new double[size];

    memset(components, 0.0, sizeof(double)*NEUMAXSQ);
//    double m_real[size][size]; double m_imag[size][size];
    double m_real[dim][dim]; double m_imag[dim][dim];
    
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            m_real[i][j] = gsl_matrix_complex_get(m,i,j).dat[0];
            m_imag[i][j] = gsl_matrix_complex_get(m,i,j).dat[1];
        }
    }

    // the following rules are valid ONLY when the initial
    // matrix is hermitian
    
    switch (dim){
        case 2:
            #include "MatrixToSU2.txt"
            break;
        case 3:
            #include "MatrixToSU3.txt"
            break;
        case 4:
            #include "MatrixToSU4.txt"
            break;
        case 5:
            #include "MatrixToSU5.txt"
            break;
        default:
            printf("GLS_MATRIX_COMPLEX to SU_vector :: Error. \n");
            exit(0);
    }
    
    //delete m_real;
    //delete m_imag;
};

SU_vector::~SU_vector(){
  //delete [] components;
};

SU_vector SU_vector::Rescale(double x){
    for(int i=0; i< SQR(dim); i++){
        components[i] = x*components[i];
    };
    return *this;
};

SU_vector SU_vector::Rotate(int i, int j, double th, double del){
    SU_vector suv = *this;
    SU_vector suv_rot(dim);
    
    switch (dim){
        case 2:
            if (i == 1 and j == 2){
                #include "RotationSU2_12.txt"
            };
            break;
        case 3:
            switch (i){
                case 1:
                    switch (j){
                        case 2:
                            #include "RotationSU3_12.txt"
                            break;
                        case 3:
                            #include "RotationSU3_13.txt"
                            break;
                    };
                    break;
                case 2:
                    switch (j){
                        case 3:
                            #include "RotationSU3_23.txt"
                            break;
                    };
            };
            break;
        case 4:
            switch (i){
                case 1:
                    switch (j){
                        case 2:
                            #include "RotationSU4_12.txt"
                            break;
                        case 3:
                            #include "RotationSU4_13.txt"
                            break;
                        case 4:
                            #include "RotationSU4_14.txt"
                            break;
                    };
                    break;
                case 2:
                    switch (j){
                        case 3:
                            #include "RotationSU4_23.txt"
                            break;
                        case 4:
                            #include "RotationSU4_24.txt"
                            break;
                    };
                    break;  
                case 3:
                    switch (j){
                        case 4:
                            #include "RotationSU4_34.txt"
                            break;
                    };
                    break;
            }
            break;
        case 5:
            switch (i){
                case 1:
                    switch (j){
                        case 2:
                            #include "RotationSU5_12.txt"
                            break;
                        case 3:
                            #include "RotationSU5_13.txt"
                            break;
                        case 4:
                            #include "RotationSU5_14.txt"
                            break;
                        case 5:
                            #include "RotationSU5_15.txt"
                            break;                            
                    };
                    break;
                case 2:
                    switch (j){
                        case 3:
                            #include "RotationSU5_23.txt"
                            break;
                        case 4:
                            #include "RotationSU5_24.txt"
                            break;
                        case 5:
                            #include "RotationSU5_25.txt"
                            break;
                    };
                    break;  
                case 3:
                    switch (j){
                        case 4:
                            #include "RotationSU5_34.txt"
                            break;
                        case 5:
                            #include "RotationSU5_35.txt"
                            break;
                    };
                    break;
                case 4:
                    switch (j){
                        case 5:
                            #include "RotationSU5_45.txt"
                            break;
                    };
                    break;
            }
            break;
        default:
            printf("SUN_rotation error. \n");
            exit(0);
    };
    
    for(int i=0; i < SQR(dim); i++){
        components[i] = suv_rot.components[i];
    };
    
    return *this;
};

SU_vector SU_vector::RotateToMass(const PhysConst* param){
    SU_vector suv = *this;
    SU_vector suv_rot(dim);
    
    double th12,th13,th23,del12,del13,del23;
    double th14,th24,th34,del14,del24,del34;
    //double th15,th25,th35,th45,del15,del25,del35,del45;
    
    switch (dim){
        case 2:
            th12 = param->th12;
            del12 = 0.0;
            #include "RotationSU2.txt"
            break;
        case 3:
            th12 = param->th12;
            del12 = 0.0;
            th13 = param->th13;
            del13 = param->delta1;
            th23 = param->th23;
            del23 = 0.0;
            #include "RotationSU3.txt"
            break;
        case 4:
            suv.Rotate(3,4,param->th34,0.0);
            suv.Rotate(2,4,param->th24,0.0);
            suv.Rotate(1,4,param->th14,param->delta2);
            suv.Rotate(2,3,param->th23,0.0);
            suv.Rotate(1,3,param->th13,param->delta1);
            suv.Rotate(1,2,param->th12,0.0);
            
            suv_rot = suv;
            
            //#include "RotationSU4.txt"
            break;
        case 5:

            break;
        default:
            printf("SUN_rotation error. \n");
            exit(0);
    };
    
    for(int i=0; i < SQR(dim); i++){
        components[i] = suv_rot.components[i];
    };
    
    return *this;
};


bool SU_vector::operator==(const SU_vector &other){
    bool equal = true;
    for(int i=0; i < SQR(dim); i++){
        if (components[i] != other.components[i]){
            equal = false;
            break;
        }
    };
    if (dim != other.dim)
        equal = false;
    
    return equal;
};

double SU_vector::operator*(const SU_vector &other){
    return SUTrace(*this,other);
};

SU_vector SU_vector::operator*(const double x){
    SU_vector su_new(dim);
    for(int i=0; i < SQR(dim); i++){
        su_new.components[i] = x*components[i];
    };
    return su_new;
};

SU_vector SU_vector::operator=(const SU_vector &other){
    dim = other.dim;
    for(int i=0; i < SQR(dim); i++){
        components[i] = other.components[i];
    };
    return *this;
};

SU_vector SU_vector::operator+=(const SU_vector &other){
    for(int i=0; i < SQR(dim); i++){
        components[i] += other.components[i];
    };
    return *this;
};

SU_vector SU_vector::operator-=(const SU_vector &other){
    for(int i=0; i < SQR(dim); i++){
        components[i] -= other.components[i];
    };
    return *this;
};

SU_vector SU_vector::operator+(const SU_vector &other){
    SU_vector su_new(dim);
    for(int i=0; i < SQR(dim); i++){
        su_new.components[i] = components[i] + other.components[i];
    };
    return su_new;
};

SU_vector SU_vector::operator-(const SU_vector &other){
    SU_vector su_new(dim);
    for(int i=0; i < SQR(dim); i++){
        su_new.components[i] = components[i] - other.components[i];
    };
    return su_new;
};

void SU_vector::InitSU_vector(gsl_matrix_complex* m){
    dim = m->size1;
    int size = SQR(dim);
    // create and set components to zero
    //components = new double[size];
    
    memset(components, 0.0, sizeof(double)*NEUMAXSQ);
//    double m_real[size][size]; double m_imag[size][size];
    double m_real[dim][dim]; double m_imag[dim][dim];
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            m_real[i][j] = gsl_matrix_complex_get(m,i,j).dat[0];
            m_imag[i][j] = gsl_matrix_complex_get(m,i,j).dat[1];
        }
    }
    // the following rules are valid ONLY when the initial
    // matrix is hermitian
    
    switch (dim){
        case 2:
            #include "MatrixToSU2.txt"
            break;
        case 3:
            #include "MatrixToSU3.txt"
            break;
        case 4:
            #include "MatrixToSU4.txt"
            break;
        case 5:
            #include "MatrixToSU5.txt"
            break;
        default:
            printf("GLS_MATRIX_COMPLEX to SU_vector :: Error. \n");
            exit(0);
    }
    //delete m_real;
    //delete m_imag;
};

void SU_vector::InitSU_vector(int d){
    dim = d;
    //components = new double[size];
    memset(components, 0.0, sizeof(double)*NEUMAXSQ);
};

/*
-----------------------------------------------------------------------
SU_vector tools
-----------------------------------------------------------------------         
*/

SU_vector operator *(const double x, SU_vector &suv){
    return suv*x;
}

SU_vector SUEvolve(SU_vector suv1, SU_vector suv2, double t){
    // here we will asumme that suv1 is the
    // evolution generator and suv2 is the evolved
    // operator. The given expressions are only valid
    // when suv1 is a projector-like SU-vector.
    int dim = suv1.dim;
    SU_vector suv_new(dim);
    switch (dim){
        case 2:
            #include "EvolutionSU2.txt"
            break;
        case 3:
            #include "EvolutionSU3.txt"
            break;
        case 4:
            #include "EvolutionSU4.txt"
            break;
        case 5:
            #include "EvolutionSU5.txt"
            break;
        default:
            printf("SUEvolution :: Error. \n");
            exit(0);
    };
    return suv_new;
}

SU_vector SUiConmutator(SU_vector suv1,SU_vector suv2){
    int dim = suv1.dim;
    SU_vector suv_new(dim);
    suv_new.dim = dim;

/*
  std::cout << "in conm " << dim << std::endl;
  std::cout << suv1 << std::endl;
  std::cout << suv2 << std::endl;
  std::cout << suv_new << std::endl;
*/
    switch (dim){
        case 2:
            #include "iConmutatorSU2.txt"
            break;
        case 3:
            #include "iConmutatorSU3.txt"
            break;
        case 4:
            #include "iConmutatorSU4.txt"
            break;
        case 5:
            #include "iConmutatorSU5.txt"
            break;
        default:
            printf("SUiConmutator :: Error. \n");
            exit(0);
    };
//  std::cout << suv_new << std::endl;
//  exit(0);
    return suv_new;
};

SU_vector SUAnticonmutator(const SU_vector suv1,const SU_vector suv2){
    int dim = suv1.dim;
    SU_vector suv_new(dim);
    switch (dim){
        case 2:
            #include "AnticonmutatorSU2.txt"
            break;
        case 3:
            #include "AnticonmutatorSU3.txt"
            break;
        case 4:
            #include "AnticonmutatorSU4.txt"
            break;
        case 5:
            #include "AnticonmutatorSU5.txt"
            break;
        default:
            printf("SUAnticonmutator :: Error. \n");
            exit(0);
    };
    return suv_new;
};

double SUTrace(SU_vector* suv1,SU_vector* suv2){
    double gen_trace = 0.0;
    double id_trace = 0.0;
    for(int i=1; i < SQR(suv1->dim); i++){
      gen_trace += (suv1->components[i])*(suv2->components[i]);
    };
    id_trace = (suv1->components[0])*(suv2->components[0])*double(suv1->dim);
    
    return id_trace+2.0*gen_trace;
};

double SUTrace(SU_vector suv1,SU_vector suv2){
    double gen_trace = 0.0;
    double id_trace = 0.0;
    
    for(int i=1; i < SQR(suv1.dim); i++){
      gen_trace += (suv1.components[i])*(suv2.components[i]);
    };
    
    id_trace = (suv1.components[0])*(suv2.components[0])*double(suv1.dim);
    return id_trace+2.0*gen_trace;
};

void SUPrint(const SU_vector* suv){
    for(int i=0; i< SQR(suv->dim); i++){
      std::cout << suv->components[i] << " ";  
    };
    std::cout << std::endl;
};

void SUPrint(const SU_vector suv){
    for(int i=0; i< SQR(suv.dim); i++){
      std::cout << suv.components[i] << " ";  
    };
    std::cout << std::endl;
};


ostream& operator<<(ostream& os, const SU_vector& V){
  int size = SQR(V.dim);
  for(int i=0; i< size-1; i++)
    os << V.components[i] << "  ";
  os <<V.components[size-1];
  return os;
}
