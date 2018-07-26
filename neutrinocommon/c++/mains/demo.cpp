/*
Author  : C.A.Arguelles
Date    : 3/JUL/2012

This files explains the basic usage of the neutrinocommon C++
Runge Kutta.

*/

#include <vector>
#include <iostream>
#include "body.h"
#include "physconst.h"
#include "neuosc.h"
#include "tools.h"

using namespace std;

int main()
{

/*
  This creates the analogus physicsconstants
  object.
*/

PhysConst* pc = new PhysConst();
/*
pc->numneu = 4;

pc->th12 = 0.563942;
pc->th13 = 0.154085;
pc->th23 = 3.14159265/4.0;

pc->dm21sq = 7.65e-5;
pc->dm31sq = 2.47e-3;

pc->delta1 = 0.0;

pc->th14 = 0.05;
pc->th24 = 0.05;
pc->th34 = 0.05;
pc->dm41sq = 0.87;
*/

pc->numneu = 3;
pc->th12 = 0.5948;
pc->th13 = 0.0;
pc->th23 = 3.14159265/4.0;

pc->delta1 = 0.0;
pc->dm21sq = 7.92e-5;
pc->dm31sq = 2.6e-3;

/*
  Lets print the current mixing parameters
*/
    cout << "Mixing angles :" << endl;
    cout << "TH12 : " << pc->th12 << endl;
    cout << "TH13 : " << pc->th13 << endl;
    cout << "TH23 : " << pc->th23 << endl;
    cout << "Square mass differences" << endl;
    cout << "DM21SQ : " << pc->dm21sq << endl;
    cout << "DM31SQ : " << pc->dm31sq << endl;
/*
  We can calculate the square mass difference
  in the flavor base corresponding to this
  parameter set.
*/

pc->Refresh();

gsl_matrix_complex* fM2;
fM2 = flavorM2(pc);
gsl_matrix_complex* U;
U = MixMatrix(pc);

printf("FM2 \n");
gsl_matrix_complex_fprintf(stdout,fM2,"%g"); 
printf("\n");

printf("U \n");
gsl_matrix_complex_fprintf(stdout,U,"%g"); 
printf("\n");

double psi[2*pc->numneu];
double EE = 50.0*(pc->GeV);

SunASnu* SA = new SunASnu;
SunASnu::Track* SA_track = new SunASnu::Track(0.0,SA->radius*pc->km/3.0);

CalNeuOscGSL(psi,pc->muon,EE,SA_track,SA,fM2,pc,1.0e-8,1.0e-8,1,0);

cout << "SunAtm State" << endl;
for(int i =0;i<2*(pc->numneu);i++){
    cout << psi[i] << endl;    
}



exit(0);


/*
  We must define the body on which the neutrino
  propagates and its corresponding trayectory.
*/

/*
 For vaccum
*/

Vacuum* vacuum = new Vacuum;
Vacuum::Track* vacuum_track = new Vacuum::Track(0.0,100.0*pc->km);

cout << "Track information" << endl;
cout << vacuum_track->x << endl;
cout << vacuum_track->xini << endl;
cout << vacuum_track->xend << endl;

cout << "Body Information" << endl;
cout << vacuum->name << endl;
cout << vacuum->density(vacuum_track) << endl;

/*
  We have to create an array where the
  calculation will be stored.
*/

double prob[pc->numneu];

/*
  Setting the neutrino energy in eV.
*/

double E = 100.0*(pc->MeV);

/*
  We calculate the oscillarion probability.
    CalNeuOscGSL(probability,
                initial_neutrino_flavor,
                neutrino_energy,
                neutrino_trayectory,
                body,
                square_mass_matrix in flavor basis,
                physics parameters,
                absolute error,
                relative error,
                return_neutrino_state,
                use optimizations)
*/
CalNeuOscGSL(prob,0,E,vacuum_track,vacuum,fM2,pc,1.0e-6,1.0e-6,0,0);
/*
  Printing the oscillation probabilities.
*/

cout << "Vacuum Oscillation Probabilities" << endl;
for(int i =0;i<pc->numneu;i++){
    cout << prob[i] << endl;    
}

/*
  We can make similar calculations for other
  bodies.
*/

/* Earth */
double baseline = 500.0*pc->km;

Earth* earth = new Earth;
Earth::Track* earth_track = new Earth::Track(0.0,baseline,baseline);

CalNeuOscGSL(prob,0,E,earth_track,earth,fM2,pc,1.0e-6,1.0e-6,0,0);

cout << "Earth Oscillation Probabilities" << endl;
for(int i =0;i<pc->numneu;i++){
    cout << prob[i] << endl;    
}

/* EarthAtm */

/* Earth */
double phi = 3.1415926;
E = 50.0*pc->GeV;

EarthAtm* earthatm = new EarthAtm;
EarthAtm::Track* earthatm_track = new EarthAtm::Track(phi);

CalNeuOscGSL(prob,1,E,earthatm_track,earthatm,fM2,pc,1.0e-6,1.0e-6,0,0);

cout << "EarthAtm Oscillation Probabilities" << endl;
for(int i = 0;i<pc->numneu;i++){
    cout << prob[i] << endl;    
}

/* Constant Density */
baseline = 10000.0*pc->km;
E = 100.0*pc->MeV;

ConstantDensity* constdens = new ConstantDensity(20.0,0.5);
ConstantDensity::Track* constdens_track = new ConstantDensity::Track(0.0,baseline);

CalNeuOscGSL(prob,0,E,constdens_track,constdens,fM2,pc,1.0e-6,1.0e-6,0,0);

cout << "Constant Density Oscillation Probabilities" << endl;
for(int i =0;i<pc->numneu;i++){
    cout << prob[i] << endl;    
}

/* Variable Density */

VariableDensity* vardens = new VariableDensity(10.0,2.0,0.5);
VariableDensity::Track* vardens_track = new VariableDensity::Track(0.0,baseline);

CalNeuOscGSL(prob,0,E,vardens_track,vardens,fM2,pc,1.0e-6,1.0e-6,0,0);

cout << "Variable Density Oscillation Probabilities" << endl;
for(int i =0;i<pc->numneu;i++){
    cout << prob[i] << endl;    
}

/* Star Density */
double star_radius = 1.0e5*pc->km;
double ini_radius = 1000.0*pc->km;
double ini_density = 4.0;

Star* star = new Star(star_radius,ini_radius,ini_density);
Star::Track* star_track = new Star::Track(0.0,star_radius);

CalNeuOscGSL(prob,0,E,star_track,star,fM2,pc,1.0e-5,0.0,0,0);

cout << "Star Oscillation Probabilities" << endl;
for(int i =0;i<pc->numneu;i++){
    cout << prob[i] << endl;    
}

/* Sun */

cout << "Sun Test" << endl;

Sun* sun = new Sun;
cout << "Getting interpolated density" << endl;
cout << sun->rdensity(0.5) << endl;
cout << "Getting interpolated hydrogen fraction" << endl;
cout << sun->rxh(0.5) << endl;

E = 100.0*pc->MeV;

Sun::Track* sun_track = new Sun::Track(0.0,sun->radius*pc->km);
CalNeuOscGSL(prob,0,E,sun_track,sun,fM2,pc,1.0e-6,1.0e-6,0,0);

cout << "Sun Oscillation Probabilities" << endl;
for(int i =0;i<pc->numneu;i++){
    cout << prob[i] << endl;    
}

/* Sun */

cout << "SunASnu Test" << endl;

SunASnu* sunatm = new SunASnu;

cout << "Getting interpolated density" << endl;
cout << sunatm->rdensity(0.75) << endl;
cout << "Getting interpolated hydrogen fraction" << endl;
cout << sunatm->rxh(0.5) << endl;

E = 10000.0*pc->MeV;

SunASnu::Track* sunatm_track = new SunASnu::Track(0.0,0.0);
CalNeuOscGSL(prob,0,E,sunatm_track,sunatm,fM2,pc,1.0e-6,1.0e-6,0,0);

cout << "SunAtm Oscillation Probabilities" << endl;
for(int i =0;i<pc->numneu;i++){
    cout << prob[i] << endl;    
}

}
