#include <vector>
#include <iostream>
#include "body.h"
#include "physconst.h"
#include "interface.h"

using namespace std;

int main()
{

/* Checking Mask */

PhysConst* pc = new PhysConst();

/*
double body_params[] = {};
double track_params[] = {0,100.0*pc.km,100.0*pc.km};
double paramlist[] = {3,
                    0.1,0.1,0.1,
                    0.0,0.0,0.0,0.0,0.0,0.0,
                    0.0,0.0,0.0,0.0,0.0,0.0,
                    1.0e-5,1.0e-3,0.0,0.0,0.0,
                    0.0,0.0,0.0};
                    
vector<double> vprob;
double testprob[pc.numneu];
                    
cout << "Mask Test" << endl;
vprob = CalNeuOscGSLMask(0,E,4,body_params,track_params,0,paramlist,1.0e-6,1.0e-6,0,0);


for(int i =0;i<pc.numneu;i++){
    //cout << testprob[i] << endl;
    cout << vprob[i] << endl;
}
*/

/*
E = 100.0*pc->MeV;

double body_params[] = {20.0,0.5};
double track_params[] = {0,10000.0*pc->km};
double paramlist[] = {3,
                    0.563942,0.154085,0.785398,
                    0.0,0.0,0.0,0.0,0.0,0.0,
                    0.0,0.0,0.0,0.0,0.0,0.0,
                    7.65e-05,0.00247,0.0,0.0,0.0,
                    0.0,0.0,0.0};
                    
vector<double> vprob;
double testprob[pc->numneu];
                    
cout << "Mask Test" << endl;
vprob = CalNeuOscGSLMask(0,E,1,body_params,track_params,0,paramlist,1.0e-6,1.0e-6,0,0);


for(int i =0;i<pc->numneu;i++){
    cout << vprob[i] << endl;
}

*/

// checking CalNeuIntGSLMask

int flux_type = 0;
double flux_params[] = {0};
double Emin = 10.0*pc->GeV;
double Emax = 100.0*pc->GeV;
int Ediv = 2;
int body_id = 1;
double body_params[] = {20.0,0.5};
double track_params[] = {0.0,100000.0*pc->km};
int neutype = 0;
double paramlist[] = {3,
                    0.563942,0.154085,0.785398,
                    0.0,0.0,0.0,0.0,0.0,0.0,
                    0.0,0.0,0.0,0.0,0.0,0.0,
                    7.65e-05,0.00247,0.0,0.0,0.0,
                    0.0,0.0,0.0};
double abs_error = 1.0e-6;
double rel_error = 1.0e-6;

int oscillation = 0;
int attenuation = 0;
int nc_inter = 0;
int cc_inter = 0;
int tau_regeneration = 0;
int muon = 0;
int electron = 0;
int basis = 0;

/*
vector< vector<double> > state = CalNeuOscIntMask(flux_type,flux_params,
                        Emin,Emax,Ediv,
                        body_id,body_params,track_params,
                        neutype,paramlist,
                        abs_error,rel_error,
                        oscillation,attenuation,nc_inter,cc_inter,
                        tau_regeneration,muon,electron,basis);

}
*/

vector< vector<double> > state = CalNeuOscSUNMask(flux_type,flux_params,
                        Emin,Emax,Ediv,
                        body_id,body_params,track_params,
                        neutype,paramlist,
                        abs_error,rel_error,
                        oscillation,attenuation,nc_inter,cc_inter,
                        tau_regeneration,muon,electron,basis);

}
