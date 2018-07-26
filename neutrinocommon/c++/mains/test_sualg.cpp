#include <iostream>
#include "neurho.h"
#include "SUNalg.h"
#include "neualg.h"
#include "physconst.h"

using namespace std;

int main()
{
    PhysConst* pc = new PhysConst();
    pc->th12 = 0.7853;
    pc->th13 = 0.1;
    pc->th23 = 0.7853;
    pc->th24 = 0.1;
    
    //flux1.RotateToMass(pc);
    
    gsl_matrix_complex* m0 = FlavorProjector(0,4);
    gsl_matrix_complex* m1 = FlavorProjector(1,4);
    gsl_matrix_complex* m2 = FlavorProjector(2,4);
    gsl_matrix_complex* m3 = FlavorProjector(3,4);
    /*
    gsl_matrix_complex_fprintf(stdout,m2,"%g");
    cout << endl ;
    gsl_matrix_complex_fprintf(stdout,m3,"%g");
    */
    SU_vector pi0(m0);
    SU_vector pi1(m1);
    SU_vector pi2(m2);
    SU_vector pi3(m3);  
    /*
    SUPrint(pi0);
    SUPrint(pi1);
    SUPrint(pi2);
    SUPrint(pi3);
    */
    SU_vector flux1 = 5.0*pi0+2.0*pi1+1.0*pi2;
    SU_vector flux2 = 3.0*pi0+1.0*pi1+5.0*pi2;
    
    SUPrint(flux1);
    flux1.RotateToMass(pc);
    SUPrint(flux1);
    
    exit(0);
    
    gsl_matrix_complex* m = FlavorProjector(0,3);
    SU_vector pi(m);


    /*    
    SUPrint(pi);
    pi.RotateToMass(pc);
    SUPrint(pi);
    //pi.Rotate(1,2,1.0,0.0);
    //SUPrint(pi);
    //pi.Rotate(1,3,2.0,0.0);
    //SUPrint(pi);
    //pi.Rotate(2,3,3.0,0.0);
    //SUPrint(pi);
    
    
    exit(0);
    */
    
    double A_comp[9] = {3.66667, 2., 3., 0., -1.5, 5., 0., 0., -2.02073};
    SU_vector A(3,A_comp);
    double B_comp[9] = {7., 3., 8., 0, -1.5, 1., 0, 0, -2.59808};
    SU_vector B(3,B_comp);
    
    double H_comp[9] = {100.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 9.0};
    cout << "b operation test" << endl;
    SUPrint(A);
    SUPrint(B);
    A -= B;
    SUPrint(A);
    
    cout << "e operation test" << endl;
    
    SU_vector H(3,H_comp);
    
    SU_vector X = 5.0*A;
    cout << X.dim << endl;
    
    SU_vector Y = A*5.0;
    cout << Y.dim << endl;
    
    cout << A*A << endl;
    
    SUPrint(SUEvolve(H,A,100.0));

    //cout << A*A << endl;    
    SUPrint(SUiConmutator(A,B));
    SUPrint(SUAnticonmutator(A,B));    
    
    exit(0);
    
    //PhysConst* pc = new PhysConst();
    //pc->th12 = 0.7853;
    //pc->th13 = 0.0;
    //pc->th23 = 0.7853;
    pc->numneu = 3;
    pc->neutype = 0;
    


    //RotateToMass(&flux1,pc);
    flux1.RotateToMass(pc);
    //RotateToMass(&flux2,pc);
    flux1.RotateToMass(pc);
    SUPrint(flux1);
    SUPrint(flux2);
    
    SUPrint(SUiConmutator(flux1,flux2));
    //SUPrint(SUAnticonmutator(x,y));
    
    return 0;
}