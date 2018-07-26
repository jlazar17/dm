#include "taumain_simpledecay.h"
#include <vector>

using namespace std;

int main(){
    double E_nu = 100.0;
    int tauola_init = 0;
    vector<double> sol = NeutrinoEnergies(E_nu,tauola_init);
}