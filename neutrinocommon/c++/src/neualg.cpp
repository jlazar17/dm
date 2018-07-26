#include "neualg.h"

// Macros
#define SQR(x)      ((x)*(x))                        // x^2
#define SQR_ABS(x)  (SQR(creal(x)) + SQR(cimag(x)))  // |x|^2
#define POW10(x)    (exp(M_LN10*(x)))                // 10^x
#define MIN(X,Y)    ( ((X) < (Y)) ? (X) : (Y) )
#define MAX(X,Y)    ( ((X) > (Y)) ? (X) : (Y) )
#define SIGN(a,b)   ( (b) > 0.0 ? (fabs(a)) : (-fabs(a)) )
#define KRONECKER(i,j)  ( (i)==(j) ? 1 : 0 )

void InitSUBox(SUBox* caja){
    int esize = caja->esize;
    int numneu = caja->numneu;
    caja->param->Refresh();
    // set branching
    caja->br_lep = 0.1736;
    // calculate energy differences
    caja->delE[0] = 0.0;
    for(int e1 = 1; e1 < caja->esize; e1++){
        caja->delE[e1] = caja->E_range[e1] - caja->E_range[e1-1];
    }

    for(int i = 0; i < numneu; i++){
        for(int ei = 0; ei < esize; ei++){
            (caja->evol_flavor_proj)[i*esize + ei].InitSU_vector(numneu);
            (caja->evol_mass_proj)[i*esize + ei].InitSU_vector(numneu);
        }
    }
    
    // initialize everything
    InitInteractions(caja);
    InitProjectors(caja);
    InitHamiltonian(caja);

    caja->numeqn = 1*(ESIZEMAX*sizeof(double) + ESIZEMAX*sizeof(double)*SQR(NEUMAX) + sizeof(int)*ESIZEMAX);
};

void InitSUBoxNoInt(SUBox* caja){
    int esize = caja->esize;
    int numneu = caja->numneu;
    caja->param->Refresh();
    // set branching
    caja->br_lep = 0.1736;
    // calculate energy differences
    caja->delE[0] = 0.0;
    for(int e1 = 1; e1 < caja->esize; e1++){
        caja->delE[e1] = caja->E_range[e1] - caja->E_range[e1-1];
    }

    for(int i = 0; i < numneu; i++){
        for(int ei = 0; ei < esize; ei++){
            (caja->evol_flavor_proj)[i*esize + ei].InitSU_vector(numneu);
            (caja->evol_mass_proj)[i*esize + ei].InitSU_vector(numneu);
        }
    }
    
    // initialize everything
    InitProjectors(caja);
    InitHamiltonian(caja);
};

void InitProjectors(SUBox* caja){    
    for(int i = 0; i < caja->numneu; i++){
        caja->mass_projectors[i].InitSU_vector(Projector(i,caja->numneu));
        caja->flavor_projectors[i].InitSU_vector(Projector(i,caja->numneu));
        //RotateToMass(&(caja->flavor_projectors[i]),caja->param);
        caja->flavor_projectors[i].RotateToMass(caja->param);
    };
    
    // initialize evol_flavor_proj
    for(int i = 0; i < caja -> numneu; i++){
        for(int ei = 0; ei < caja -> esize; ei++){
            (caja->evol_flavor_proj)[i*caja->esize + ei] = caja->flavor_projectors[i];
            (caja->evol_mass_proj)[i*caja->esize + ei] = caja->mass_projectors[i];
        }
    }
};

void InitHamiltonian(SUBox* caja){
    caja->DM2.InitSU_vector(caja->numneu);
    for(int i = 1; i < caja->numneu; i++){
        caja->DM2 += (caja->mass_projectors[i])*gsl_matrix_get(caja->param->dmsq,i,0);
    };
};

void InitInteractions(SUBox* caja){
    int esize = caja->esize;
    int neutype = caja->neutype;
    //units
    double cm2GeV = pow(caja->param->cm,2)*pow(caja->param->GeV,-1);
    double cm2 = pow(caja->param->cm,2);
    double GeVm1 = pow(caja->param->GeV,-1);
    // load cross sections
    
    NeutrinoCrossSections ncs(caja->E_range[0],caja->E_range[esize-1],esize-1);
    // initializing cross section arrays
    double dsignudE_NC[NEUINT][ESIZEMAX][ESIZEMAX];
    double dsignudE_CC[NEUINT][ESIZEMAX][ESIZEMAX];
    
    // filling cross section arrays
    for(int flv = 0; flv < NEUINT; flv++){
        for(int e1 = 0; e1 < esize; e1++){
            // differential cross sections
            for(int e2 = 0; e2 < e1; e2++){
                dsignudE_NC[flv][e1][e2] = ncs.dsde_NC(e1,e2,flv,neutype)*cm2GeV;
                dsignudE_CC[flv][e1][e2] = ncs.dsde_CC(e1,e2,flv,neutype)*cm2GeV;
            }
            // total cross sections
            caja->sigma_CC[flv][e1] = ncs.sigma_CC(e1,flv,neutype)*cm2;
            caja->sigma_NC[flv][e1] = ncs.sigma_NC(e1,flv,neutype)*cm2;
        }
    }
    
    #ifdef FixCrossSections
    // fix charge current and neutral current differential cross sections
    double XCC_MIN,XNC_MIN,XCC_int,XNC_int,CC_rescale,NC_rescale;
    for(int flv = 0; flv < NEUINT; flv++){
        XCC_MIN = caja->sigma_CC[flv][0];
        XNC_MIN = caja->sigma_CC[flv][0];
        for(int e1 = 0; e1 < esize; e1++){
            XCC_int = 0.0;
            XNC_int = 0.0;
            for(int e2 = 0; e2 < e1; e2++){
                XCC_int += dsignudE_CC[flv][e1][e2]*caja->delE[e2];
                XNC_int += dsignudE_NC[flv][e1][e2]*caja->delE[e2];
            }
            
            if(e1 != 0 ){
                CC_rescale = (caja->sigma_CC[flv][e1] - XCC_MIN)/XCC_int;
                NC_rescale = (caja->sigma_NC[flv][e1] - XNC_MIN)/XNC_int;
                
                for(int e2 = 0; e2 < e1; e2++){
                    dsignudE_CC[flv][e1][e2] = dsignudE_CC[flv][e1][e2]*CC_rescale;
                    dsignudE_NC[flv][e1][e2] = dsignudE_NC[flv][e1][e2]*NC_rescale;
                }
            }
        }
    }
    #endif

    // constructing dNdE    
    for(int flv = 0; flv < NEUINT; flv++){
        for(int e1 = 0; e1 < esize; e1++){
            for(int e2 = 0; e2 < e1; e2++){
                if (dsignudE_NC[flv][e1][e2] < 1.0e-50){
                    caja->dNdE_NC[flv][e1][e2] = 0.0;
                } else {
                    caja->dNdE_NC[flv][e1][e2] = (dsignudE_NC[flv][e1][e2])/(caja->sigma_NC[flv][e1]);
                }
                if (dsignudE_CC[flv][e1][e2] < 1.0e-50){
                    caja->dNdE_CC[flv][e1][e2] = 0.0;
                } else {
                    caja->dNdE_CC[flv][e1][e2] = (dsignudE_CC[flv][e1][e2])/(caja->sigma_CC[flv][e1]);
                }           
            }
        }
    }

    // initialize interaction lenghts to zero    
    // tau decay length array
    for(int e1 = 0; e1 < esize; e1++){
        caja->invlen_tau[e1] = 1.0/((caja->param->tau_lifetime)*(caja->E_range[e1]*caja->param->tau_mass));
    }
    
    // load tau decay spectra
    TauDecaySpectra tdc(caja->E_range[0],caja->E_range[esize-1],esize-1);
    // fill arrays
    double dNdENuTau_All[ESIZEMAX][ESIZEMAX];
    double dNdENuTau_Lep[ESIZEMAX][ESIZEMAX];
    double dNdELeTau_All[ESIZEMAX][ESIZEMAX];
    double dNdELeTau_Lep[ESIZEMAX][ESIZEMAX];
    
    for(int e1 = 0; e1 < esize; e1++){
        for(int e2 = 0; e2 < e1; e2++){
            dNdENuTau_All[e1][e2] = tdc.dNdEnu_All(e1,e2)*GeVm1;
            dNdENuTau_Lep[e1][e2] = tdc.dNdEnu_Lep(e1,e2)*GeVm1;
            dNdELeTau_All[e1][e2] = tdc.dNdEle_All(e1,e2)*GeVm1;
            dNdELeTau_Lep[e1][e2] = tdc.dNdEle_Lep(e1,e2)*GeVm1;
        }
    }
    
    // constructing dNdE_tau_lep
    for(int e1 = 0; e1 < esize; e1++){
        for(int e2 = 0; e2 < e1; e2++){
            caja->dNdE_tau_all[e1][e2] = dNdENuTau_All[e1][e2];
            caja->dNdE_tau_lep[e1][e2] = dNdENuTau_Lep[e1][e2];
        }
    }
    
    #ifdef FixCrossSections
    // fix tau decay spectra cross section
    double tau_all_int, tau_lep_int,tau_lep_rescale,tau_all_rescale;
    for(int e1 = 0; e1 < esize; e1++){
        tau_all_int = 0.0;
        tau_lep_int = 0.0;
        for(int e2 = 0; e2 < e1; e2++){
             tau_all_int += caja->dNdE_tau_all[e1][e2]*caja->delE[e2];
             tau_lep_int += caja->dNdE_tau_lep[e1][e2]*caja->delE[e2];
        }
        
        if( caja->dNdE_tau_all[e1][0]*caja->E_range[0] < 0.25 ) {
            tau_all_rescale = (1.0 - caja->dNdE_tau_all[e1][0]*caja->E_range[0])/tau_all_int;
            tau_lep_rescale = (caja->br_lep - caja->dNdE_tau_lep[e1][0]*caja->E_range[0])/tau_lep_int;
            
            for(int e2 = 0; e2 < e1; e2++){
                caja->dNdE_tau_all[e1][e2] = caja->dNdE_tau_all[e1][e2]*tau_all_rescale;
                caja->dNdE_tau_lep[e1][e2] = caja->dNdE_tau_lep[e1][e2]*tau_lep_rescale;
            }
        }
    }
    #endif
};

void UpdateInteractions(SUBox* caja){
    double density = caja->body->density(caja->track);
    //double num_nuc = caja->param->Na*pow(caja->param->cm,-3)*density;
    double num_nuc = (caja->param->gr*pow(caja->param->cm,-3))*density*2.0/(caja->param->proton_mass+caja->param->neutron_mass);
    
    #ifdef UpdateInteractions_DEBUG
        cout << "Density " << density << endl;
        cout << "Nucleon Number " << num_nuc << endl;
    #endif
    
    if(num_nuc == 0){
        num_nuc = caja->param->Na*pow(caja->param->cm,-3)*1.0e-10;
    }
    
    for(int flv = 0; flv < NEUINT; flv++){
            #ifdef UpdateInteractions_DEBUG
            cout << "============" << flv << "============" << endl;
            #endif
        for(int e1 = 0; e1 < caja->esize; e1++){
            #ifdef UpdateInteractions_DEBUG
                cout << "== CC NC Terms x = " << caja->track->x/caja->param->km << " [km] ";
                cout << "E = " << caja->E_range[e1] << " [eV] ==" << endl;
                cout << "CC : " << caja->sigma_CC[flv][e1]*num_nuc << " NC : " << caja->sigma_NC[flv][e1]*num_nuc << endl;
                cout << "==" << endl;
            #endif
            caja->invlen_NC[flv][e1] = caja->sigma_NC[flv][e1]*num_nuc;
            caja->invlen_CC[flv][e1] = caja->sigma_CC[flv][e1]*num_nuc;
            caja->invlen_int[flv][e1] = caja->invlen_NC[flv][e1] + caja->invlen_CC[flv][e1];
        }
    }
};

/*
void RotateToFlavor(SU_vector* suv,PhysConst* param){
    // performs inverse rotation
    switch (suv->dim){
        case 2:
            suv->Rotate(1,2,-param->th12,0.0);
            break;
        case 3:
            suv->Rotate(2,3,-param->th23,0.0);
            suv->Rotate(1,3,-param->th13,-param->delta1);
            suv->Rotate(1,2,-param->th12,0.0);
            break;
        case 4:
            suv->Rotate(3,4,-param->th34,0.0);
            suv->Rotate(2,4,-param->th24,0.0);
            suv->Rotate(1,4,-param->th14,-param->delta2);
            suv->Rotate(2,3,-param->th23,0.0);
            suv->Rotate(1,3,-param->th13,-param->delta1);
            suv->Rotate(1,2,-param->th12,0.0);
            break;
        case 5:
            suv->Rotate(4,5,-param->th45,0.0);
            suv->Rotate(3,5,-param->th35,0.0);
            suv->Rotate(2,5,-param->th25,0.0);
            suv->Rotate(1,5,-param->th15,0.0);
            suv->Rotate(3,4,-param->th34,0.0);
            suv->Rotate(2,4,-param->th24,0.0);
            suv->Rotate(1,4,-param->th14,-param->delta2);
            suv->Rotate(2,3,-param->th23,0.0);
            suv->Rotate(1,3,-param->th13,-param->delta1);
            suv->Rotate(1,2,-param->th12,0.0);
            break;
    };
};
*/

void EvolveProjectors(double x,SUBox* caja){
    int esize = caja->esize;
    
    for(int i = 0; i < caja -> numneu; i++){
        for(int ei = 0; ei < caja -> esize; ei++){
            (caja->evol_flavor_proj)[i*esize + ei] = SUEvolve((caja->DM2)*(0.5/caja->E_range[ei]),caja->flavor_projectors[i],x);
            (caja->evol_mass_proj)[i*esize + ei] = SUEvolve((caja->DM2)*(0.5/caja->E_range[ei]),caja->mass_projectors[i],x);
        }
    }
};

SU_vector SUHInt(int ei,SUBox* caja){

    SU_vector* projectors = caja->evol_flavor_proj;
    int esize = caja->esize;
    double ye = caja->body->ye(caja->track);
    double density = caja->body->density(caja->track);


    double CC = caja->param->sqrt2*caja->param->GF*caja->param->Na*pow(caja->param->cm,-3)*density*ye;
    double NC = CC*(-0.5*(1.0-ye)/ye);

    /*
    std::cout << "density" << std::endl;
    std::cout << density << " " << ye << " " << caja->param->sqrt2*caja->param->GF*caja->param->Na*pow(caja->param->cm,-3) << " " << CC <<" "<< NC << std::endl;
    std::cout << caja->param->sqrt2 << " " << caja->param->GF << " " << caja->param-> Na << " " << pow(caja->param->cm,-3) << endl;
    */
    // build SU_vector
    SU_vector potential(caja->numneu);
    // construct potential in flavor basis
    potential = (CC+NC)*projectors[0*esize + ei] + (NC)*projectors[1*esize + ei] + (NC)*projectors[2*esize + ei];

    //std::cout << caja->track->x/caja->param->km << " ";
    //std::cout << potential << std::endl;

    if (caja->neutype == 0){
        return potential;
    } else if (caja->neutype == 1){
        return (-1.0)*potential;
    } else{
        exit(0);
    }
};

void TauLeptonReinj_SUN(double *state_dbl_neu, void* box_neu){
    // format container
    SUBox* caja_neu = (SUBox*) box_neu;
    // state vector box
    SUState* su_state_neu = (SUState*) state_dbl_neu;
    // vector format
    SU_vector* suv_array_neu = su_state_neu->neu_array;
    double* lep_array_neu = su_state_neu->lep_array;
    // branching
    int esize = caja_neu->esize;
    
    for(int e1 = 0; e1 < esize; e1++){
        double tau_neu_all  = 0.0;
        
        for(int e2 = e1 +1; e2 < esize; e2++){
            tau_neu_all += caja_neu->dNdE_tau_all[e2][e1]*caja_neu->delE[e2]*lep_array_neu[e2];
        }
        // adding new fluxes   
        suv_array_neu[e1] += tau_neu_all*caja_neu->evol_flavor_proj[2*esize+e1];
    }
    // clean all lepton arrays
    for(int e1 = 0; e1 < esize; e1++){
        lep_array_neu[e1] = 0.0;
    }
};

void TauLeptonReinj_SUN(double *state_dbl_neu, void* box_neu, double * state_dbl_aneu, void *box_aneu){
    // format container
    SUBox* caja_neu = (SUBox*) box_neu;
    SUBox* caja_aneu = (SUBox*) box_aneu;
    // state vector box
    SUState* su_state_neu = (SUState*) state_dbl_neu;
    SUState* su_state_aneu = (SUState*) state_dbl_aneu;
    // vector format
    SU_vector* suv_array_neu = su_state_neu->neu_array;
    SU_vector* suv_array_aneu = su_state_aneu->neu_array;
    double* lep_array_neu = su_state_neu->lep_array;
    double* lep_array_aneu = su_state_aneu->lep_array;
    // branching
    int esize = caja_neu->esize;
    // cross section cosos
    
    for(int e1 = 0; e1 < esize; e1++){
        double tau_neu_all  = 0.0;
        double tau_neu_lep  = 0.0;
        double tau_aneu_all = 0.0;
        double tau_aneu_lep = 0.0;
        
        for(int e2 = e1 +1; e2 < esize; e2++){
        #ifdef TauLeptonReinj_SUN_DEBUG
            cout << "TauNeuAll DifXsections" << endl;
            cout << e1 << " " << e2 << " "<< caja_neu->dNdE_tau_all[e2][e1] << endl;
            cout << "TauFlux" << endl;
            cout << e2 << " "<< lep_array_neu[e2] << endl;
        #endif
            tau_neu_all  += caja_neu->dNdE_tau_all[e2][e1]*caja_neu->delE[e2]*lep_array_neu[e2];
            tau_neu_lep  += caja_neu->dNdE_tau_lep[e2][e1]*caja_neu->delE[e2]*lep_array_neu[e2];
            tau_aneu_all += caja_aneu->dNdE_tau_all[e2][e1]*caja_aneu->delE[e2]*lep_array_aneu[e2];
            tau_aneu_lep += caja_aneu->dNdE_tau_lep[e2][e1]*caja_aneu->delE[e2]*lep_array_aneu[e2];
        }
        // note that the br_lepton is already included in dNdE_tau_lep
        // adding new fluxes
        suv_array_neu[e1]  += tau_neu_all*caja_neu->evol_flavor_proj[2*esize+e1] +
                              tau_aneu_lep*(caja_neu->evol_flavor_proj[0*esize+e1])+
                              tau_aneu_lep*(caja_neu->evol_flavor_proj[1*esize+e1]);
        suv_array_aneu[e1] += tau_aneu_all*caja_aneu->evol_flavor_proj[2*esize+e1] +
                              tau_neu_lep*(caja_aneu->evol_flavor_proj[0*esize+e1])+
                              tau_neu_lep*(caja_aneu->evol_flavor_proj[1*esize+e1]);
    }
    
    // clean all lepton arrays

    for(int e1 = 0; e1 < esize; e1++){
        lep_array_neu[e1] = 0.0;
        lep_array_aneu[e1] = 0.0;
    }

};

int RHS_SUN(double x ,const double *state_dbl_in,double *state_dbl_out,void *box){
    // format container
    SUBox* caja = (SUBox*) box;
    int esize = caja->esize;
    // format state structure
    SUState* state_in = (SUState*) state_dbl_in;
    SU_vector* suv_array_in = state_in->neu_array;
    //double* lep_array_in = state_in->lep_array;
    // create output structure
    SUState* state_out = (SUState*) state_dbl_out;
    SU_vector* suv_array_out = state_out->neu_array;
    double* lep_array_out = state_out->lep_array;
    
    #ifdef RHS_SUN_POS_DEBUG
        cout << "x : " << x << endl;
    #endif
    
    // change location
    caja->track->x = x;
    if(caja->oscillation){
        EvolveProjectors(x,caja);  
    }
    
    // update interactions
    if(caja->nc_inter or caja->attenuation){
        UpdateInteractions(caja);
    }

    //double br_lep = caja->br_lep;
    for(int e1 = 0; e1 < esize; e1++){
        // fix mierda
        suv_array_in[e1].dim = caja->numneu;
        suv_array_out[e1] = SU_vector(caja->numneu);
    }
    
    #ifdef RHS_SUN_NC_EXT_DEBUG
        SU_vector nc_term_sub(caja->numneu);
        SU_vector nc_term_add(caja->numneu);
    #endif
   /* 
    cout.precision(15);
    cout << " x " << x << " rho " << suv_array_in[0] << endl;
    cout << " x " << x << " drho ";
    SUPrint(SUiConmutator(suv_array_in[0],SUHInt(0,caja)));
    cout << " x " << x << " hi " << SUHInt(0,caja) << endl;
*/
/*
    exit(1);
*/
    // neutrino part
    for(int e1 = 0; e1 < esize; e1++){
        // conmutator term
        if(caja->oscillation){
            #ifdef RHS_SUN_OSC_DEBUG
                cout << "SUiConmutator term" << endl;
                SUPrint(SUiConmutator(suv_array_in[e1],SUHInt(e1,caja)));
            #endif
            suv_array_out[e1] = SUiConmutator(suv_array_in[e1],SUHInt(e1,caja));
        }
        
        // anticonmutator term
        if(caja->attenuation){
            #ifdef RHS_SUN_ATT_DEBUG
                cout << "Interaction Lengths" << endl;
                cout << caja->invlen_int[0][e1] << " " << caja->invlen_int[1][e1] << " " << caja->invlen_int[2][e1] << endl;
                
                cout << "SUAnticonmutator term" << endl;
                SUPrint(SUAnticonmutator(suv_array_in[e1], caja->evol_flavor_proj[0*esize+e1]*(0.5*caja->invlen_int[0][e1]) +
                                                                        caja->evol_flavor_proj[1*esize+e1]*(0.5*caja->invlen_int[1][e1]) +
                                                                        caja->evol_flavor_proj[2*esize+e1]*(0.5*caja->invlen_int[2][e1])));
            #endif
            
            #ifdef RHS_SUN_OLD_ATT
                suv_array_out[e1] -= SUAnticonmutator(caja->evol_flavor_proj[0*esize+e1]*(0.5*caja->invlen_int[0][e1]) +
                                                      caja->evol_flavor_proj[1*esize+e1]*(0.5*caja->invlen_int[1][e1]) +
                                                      caja->evol_flavor_proj[2*esize+e1]*(0.5*caja->invlen_int[2][e1]),
                                                      suv_array_in[e1]);
            
            #else
                suv_array_out[e1] -= SUAnticonmutator(caja->evol_flavor_proj[0*esize+e1]*(0.5*caja->invlen_CC[0][e1]) +
                                                      caja->evol_flavor_proj[1*esize+e1]*(0.5*caja->invlen_CC[1][e1]) +
                                                      caja->evol_flavor_proj[2*esize+e1]*(0.5*caja->invlen_CC[2][e1]),
                                                      suv_array_in[e1]);
            #endif
        }
        
        #ifdef RHS_SUN_NC_EXT_DEBUG
            nc_term_sub += SUAnticonmutator(suv_array_in[e1], caja->evol_flavor_proj[0*esize+e1]*(0.5*caja->invlen_NC[0][e1]) +
                                                                        caja->evol_flavor_proj[1*esize+e1]*(0.5*caja->invlen_NC[1][e1]) +
                                                                        caja->evol_flavor_proj[2*esize+e1]*(0.5*caja->invlen_NC[2][e1]));
            
            cout << "interactions webas" << endl;
            cout << caja->invlen_CC[0][e1] << " " << caja->invlen_NC[0][e1] << " " << caja->invlen_int[0][e1] << endl;
        #endif

        // nc term
        if(caja->nc_inter){
            for(int e2 = e1 + 1; e2 < esize; e2++){
                #ifdef RHS_SUN_NC_DEBUG
                    cout << "NC term " << e1 << " " << e2 << endl;             
                    SUPrint(suv_array_in[e2]*((caja->dNdE_NC[0][e2][e1])*(caja->invlen_NC[0][e2])*(caja->delE[e2])));
                #endif
                
                #ifdef RHS_SUN_NC_EXT_DEBUG
                    nc_term_add += suv_array_in[e2]*((caja->dNdE_NC[0][e2][e1])*(caja->invlen_NC[0][e2])*(caja->delE[e2]));
                #endif
                
                #ifdef SterileNC
                    suv_array_out[e1] += SUAnticonmutator(caja->evol_flavor_proj[0*esize+e1] + caja->evol_flavor_proj[1*esize+e1] + caja->evol_flavor_proj[2*esize+e1],
                                                          suv_array_in[e2])*(0.5*(caja->dNdE_NC[0][e2][e1])*(caja->invlen_NC[0][e2])*(caja->delE[e2]));
                #else
                    suv_array_out[e1] += suv_array_in[e2]*((caja->dNdE_NC[0][e2][e1])*(caja->invlen_NC[0][e2])*(caja->delE[e2]));
                #endif
            }
            #ifndef RHS_SUN_OLD_ATT
                for(int e2 = 0; e2 < e1; e2 ++){
                    suv_array_out[e1] -= suv_array_in[e1]*((caja->dNdE_NC[0][e1][e2])*(caja->invlen_NC[0][e2])*(caja->delE[e2]));
                }
            #endif
        }
        
        #ifdef RHS_SUN_NC_EXT_DEBUG
            cout << "SUstate NC terms" << endl;
            cout << nc_term_add.components[2] << endl;
            SUPrint(nc_term_add);
            SUPrint(nc_term_sub);
            
            cout << nc_term_add*nc_term_add << endl;
            cout << nc_term_sub*nc_term_sub << endl;
        #endif

        /*
        // tau regeneration terms will be implemented manually
        for(int e2 = e1 +1; e2 < caja->esize; e2++){
            suv_array_out[e1] += br_lep*caja->evol_flavor_proj[2][e1]*suv_array_out[e2]*caja->dNdE_tau_lep[e2][e1]*caja->invlen_tau[e2]*caja->delE[e2];
        };
        */
        
        #ifdef RHS_SUN_OUT_DEBUG
            cout << "SUstate out" << endl;
            SUPrint(suv_array_out[e1]);
        #endif    
    }
    
    // lepton part
    if(caja->tau_regeneration){
        for(int e1 = 0; e1 < caja->esize; e1++){
            lep_array_out[e1] = 0.0;
            for(int e2 = e1 +1; e2 < caja->esize; e2++){
                #ifdef RHS_SUN_TAU_DEBUG
                    cout << "Tau term " << e1 << " " << e2 << endl;
                    cout << caja->invlen_CC[2][e2] << endl;
                    cout << caja->dNdE_CC[2][e2][e1] << endl;
                    cout << (caja->evol_flavor_proj[2*esize+e2]*suv_array_in[e2])*caja->invlen_CC[2][e2]*caja->dNdE_CC[2][e2][e1]*caja->delE[e2] << endl;
                #endif
                
                lep_array_out[e1] += (caja->evol_flavor_proj[2*esize+e2]*suv_array_in[e2])*caja->invlen_CC[2][e2]*caja->dNdE_CC[2][e2][e1]*caja->delE[e2];
            }
            
            /*
            // interaction length will be disregarted and will be make manually
                lep_array_out[e1] -= lep_array_in[e1]*invlen_tau[e1]
            */
        }     
    }
    
    state_dbl_out = (double*) state_out;
    return GSL_SUCCESS;
};


int JAC_SUN(double,const double *, double * , double [], void *){
    return GSL_SUCCESS;
};

void ConvertOscSetupToSUBox(OscSetup* setup,SUBox* caja){
    // just copy everything need
    caja->numneu = setup->param->numneu;
    caja->neutype = setup->neutype;
    caja->track = setup->track;
    caja->body = setup->body;
    caja->param = setup->param;
    caja->oscillation = setup->oscillation;
    caja->attenuation = setup->attenuation;
    caja->nc_inter = setup->nc_inter;
    caja->cc_inter = setup->cc_inter;
    caja->tau_regeneration = setup->tau_regeneration;
    
    caja->esize = (const int) setup->E_range.size();
    copy(setup->E_range.begin(),setup->E_range.end(),caja->E_range);
};

SUState* ConvertSUVectorToSUState(SUVector* sv,SUBox* caja){
    int esize = sv->size();
    // create structures
    SUState* su_state = new SUState;
    //su_state->neu_array = (SU_vector*) malloc(esize*sizeof(SU_vector));
    //su_state->lep_array = (double*) malloc(esize*sizeof(double));

    for(int e = 0; e < esize; e++){
        // neutrino part
        su_state->neu_array[e] = sv->at(e);
        // go to mass basis
        su_state->neu_array[e].RotateToMass(caja->param);
        // lepton part
        su_state->lep_array[e] = 0.0;
    }
    std::cout << "hola" << ESIZEMAX << " " << NEUMAX << std::endl;
    caja->numeqn = 1*(ESIZEMAX*sizeof(double) + ESIZEMAX*sizeof(double)*SQR(NEUMAX) + sizeof(int)*ESIZEMAX);
    
    return su_state;
};

SUState* ConvertVectorStateToSUState(VectorFluxState* vfs,SUBox* caja){
    int esize = vfs->size();
    // create structures
    SUState* su_state = new SUState;
    //su_state->neu_array = (SU_vector*) malloc(esize*sizeof(SU_vector));
    //su_state->lep_array = (double*) malloc(esize*sizeof(double));

    for(int e = 0; e < esize; e++){
        FluxState* state = vfs->at(e);
        // neutrino part
        su_state->neu_array[e].InitSU_vector(state->F_nu);
        // go to mass basis
        su_state->neu_array[e].RotateToMass(caja->param);
        // lepton part
        su_state->lep_array[e] = 0.0;
        
        //gsl_matrix_complex_fprintf(stdout,state->F_nu,"%g");
        //SUPrint(su_state->neu_array[e]);
        //exit(0);
        //free(state);
    }
    
    //caja->numeqn = 2*(esize*sizeof(double) + esize*sizeof(double)*SQR(caja->numneu) + 2*sizeof(int)*esize);
    caja->numeqn = 1*(ESIZEMAX*sizeof(double) + ESIZEMAX*sizeof(double)*SQR(NEUMAX) + sizeof(int)*ESIZEMAX);
    
    return su_state;    
};

VectorFluxState* ConvertSUStateToVectorState(SUState* su_state,SUBox* caja){
    int e_size = caja->esize;
    int numneu = caja->numneu;
    
    VectorFluxState* state_vector = new VectorFluxState;
    for(int e = 0; e < e_size; e++){
        FluxState* state = new FluxState;
        // nice way
        //state->F_nu = FromSUVectorToGSLMatrix(su_state->neu_array[e]);
        // initialize neutrino matrix
        state->F_nu = gsl_matrix_complex_calloc(numneu,numneu);
        state->FluxFlavor = new double[numneu];
        state->FluxMass = new double[numneu];
        
        for(int flv = 0; flv < numneu; flv++){
            gsl_matrix_complex_set(state->F_nu,flv,flv,
                                   gsl_complex_rect((caja->evol_flavor_proj[flv*e_size+e])*su_state->neu_array[e],
                                                    0.0)
                                   );
            // initialize real flux array
            // flavor fluxes
            state->FluxFlavor[flv] = (caja->evol_flavor_proj[flv*e_size+e])*su_state->neu_array[e];
            // mass fluxes
            state->FluxMass[flv] = (caja->evol_mass_proj[flv*e_size+e])*su_state->neu_array[e];
            // flux state
            state->FluxState = su_state->neu_array[e];
        }
        // initialize lepton matrix
        state->F_le = gsl_matrix_complex_calloc(numneu,1);
        gsl_matrix_complex_set(state->F_le,2,0,gsl_complex_rect(su_state->lep_array[e],0.0));

        //push back
        state_vector->push_back(state);
    }
    return state_vector;
};

int CalNeuOscSUN(SUVector* state_vector_in[2], SUBox* caja[2],double abs_error,double rel_error){
    int gsl_status = GSL_SUCCESS;
    bool neu_and_aneu = false;
    
    #ifdef CalNeuOscSUN_DEBUG
        printf("Configuring boxes.\n");
    #endif
    
    if(caja[0]->neutype == 0 or caja[0]->neutype == 1) {
        // one and only one neutrino type
        caja[0]->param->neutype = caja[0]->neutype;
        neu_and_aneu = false;
    } else {
        // neutrinos and antineutrinos
        neu_and_aneu = true;
    }

    // initial box setup
    if(neu_and_aneu){
        caja[0]->neutype = 0;
        caja[0]->param->neutype = 0;
        
        caja[1]->neutype = 1;
        caja[1]->param->neutype = 1;
    }

    InitSUBox(caja[0]);
    InitSUBox(caja[1]);
    
    #ifdef CalNeuOscSUN_DEBUG
        printf("Begin convert SUVector to SUState.\n");
    #endif
    
    // setup initial state    
    SUState* state_in = ConvertSUVectorToSUState(state_vector_in[0],caja[0]);

    SUState* state_in_aux = NULL;
    if (neu_and_aneu){
        state_in_aux = ConvertSUVectorToSUState(state_vector_in[1],caja[1]);
    }
    
    #ifdef CalNeuOscSUN_DEBUG
        printf("End convert SUVector to SUState.\n");
    #endif
    
    #ifdef CalNeuOscSUN_DEBUG
        printf("Convert SUState to double.\n");
    #endif

    // convert to double
    double* state_dbl = (double*) state_in;
    double* state_dbl_aux = (double*) state_in_aux;
    
    // setting up GSL ODE solver
    gsl_odeiv2_step *s = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45, caja[0]->numeqn);
    gsl_odeiv2_control *c = gsl_odeiv2_control_y_new(abs_error,rel_error);
    gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc(caja[0]->numeqn);
    
    gsl_odeiv2_step *s_aux = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45, caja[1]->numeqn);
    gsl_odeiv2_control *c_aux = gsl_odeiv2_control_y_new(abs_error,rel_error);
    gsl_odeiv2_evolve *e_aux = gsl_odeiv2_evolve_alloc(caja[1]->numeqn);
    
    // ODE system
    gsl_odeiv2_system sys = {&RHS_SUN, NULL, caja[0]->numeqn, caja[0]};
    gsl_odeiv2_system sys_aux = {&RHS_SUN, NULL, caja[1]->numeqn, caja[1]};
    
    // defining ODE extra variables
    double km = caja[0]->param->km;
    
    double x = 0;                       // ODE independent variable
    double x_aux = 0;                   // ODE2 independent variable
    double x_ini = caja[0]->track->xini;  // initial position
    double x_end = caja[0]->track->xend;  // final position
    // step sizes
    double h        = MIN(10.0*km,x_end/10.0);
    double h_min    = 1.0e-5*km;
    double h_max    = MIN(100.0*km,x_end/5.0);
    
    #ifdef USE_GSL_ODE_DRIVER
        // primary driver
        gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rkf45,abs_error,rel_error,0.0);
        gsl_odeiv2_step_set_driver(s,d);
        //gsl_odeiv2_driver_set_hmin(d,h_min);
        gsl_odeiv2_driver_set_hmax(d,h_max);
        
        // auxiliary driver
        gsl_odeiv2_driver *d_aux = gsl_odeiv2_driver_alloc_y_new(&sys_aux,gsl_odeiv2_step_rkf45,abs_error,rel_error,0.0);
        gsl_odeiv2_step_set_driver(s_aux,d_aux);
        //gsl_odeiv2_driver_set_hmin(d_aux,h_min);
        gsl_odeiv2_driver_set_hmax(d_aux,h_max);
    #endif  
    
    #ifdef CalNeuOscSUN_DEBUG
        printf("GSL paramers :\n");
        printf("x_ini : %g [km]\n", x_ini/km);
        printf("x_end : %g [km]\n", x_end/km);
        printf("h : %g [km]\n", h/km);
        printf("h_min : %g [km]\n", h_min/km);      
    #endif
    
    // initial position
    x = x_ini;
    x_aux = x_ini;
    
    #ifdef CalNeuOscSUN_DEBUG
        int count = 0;
        int count_step = 10;
    #endif
    
    #ifdef USE_FIX_STEP
        h = 1.0*km;
    #endif

    #ifdef CalNeuOscSUN_DEBUG
        printf("Start calculation.\n");
    #endif
    
    while (x < x_end){
        #ifdef USE_GSL_ODE_DRIVER
            double x_inter = x + 50.0*km;
            gsl_status = gsl_odeiv2_driver_apply(d,&x,x_inter,state_dbl);
            if( neu_and_aneu ){
                gsl_status = gsl_odeiv2_driver_apply(d_aux,&x_aux,x_inter,state_dbl_aux);
            }
        #else
            gsl_status = gsl_odeiv2_evolve_apply(e,c,s,&sys,&x,x_end,&h,state_dbl);
            if( neu_and_aneu ){
                gsl_status = gsl_odeiv2_evolve_apply(e_aux,c_aux,s_aux,&sys_aux,&x_aux,x_end,&h,state_dbl_aux);
            }
        #endif


        #ifdef ManualTauReinjection
            if(caja[0]->tau_regeneration and not neu_and_aneu){
                TauLeptonReinj_SUN(state_dbl,caja[0]);
            } else if (caja[0]->tau_regeneration)  {
                TauLeptonReinj_SUN(state_dbl,caja[0],state_dbl_aux,caja[1]);
            }
        #endif
        
        #ifdef CalNeuOscSUN_DEBUG
            if(count%count_step == 0){
                printf("x_current : %g %g [km]\n", x/km,x_aux/km);
                //printf("step : %g [km]\n", h/km);
            }
        #endif

        if( gsl_status != GSL_SUCCESS ){
            fprintf(stderr,"CalNeuOscGSL: Error in GSL ODE solver.\n");
            break;
        }
        
        if(h < h_min){
            h = h_min;
        }
        
        if(h > h_max){
            h = h_max;
        }
        
        #ifdef CalNeuOscSUN_DEBUG
            count++;
        #endif
    }
    
    #ifdef CalNeuOscSUN_DEBUG
        printf("End calculation.\n");
    #endif
    
    if (e) { gsl_odeiv2_evolve_free(e);      e = NULL; }
    if (c) { gsl_odeiv2_control_free(c);     c = NULL; }
    if (s) { gsl_odeiv2_step_free(s);        s = NULL; }
    #ifdef USE_GSL_ODE_DRIVER
    if (d) { gsl_odeiv2_driver_free(d);        d = NULL; }
    if (d_aux) { gsl_odeiv2_driver_free(d_aux);        d_aux = NULL; }    
    #endif

    // convert back to SUState
    SUState* state_out = (SUState*) state_dbl;
    SUState* state_out_aux = (SUState*) state_dbl_aux;
    
    for(int i = 0; i < caja[0]->esize; i++)
    {
        state_vector_in[0]->at(i) = state_out->neu_array[i];
        if(neu_and_aneu){
            state_vector_in[1]->at(i) = state_out_aux->neu_array[i];
        }
    }

    return GSL_SUCCESS;
}

int CalNeuOscSUN(VectorFluxState* state_vector_in[2], SUBox* caja[2],double abs_error,double rel_error){
    int gsl_status = GSL_SUCCESS;
    bool neu_and_aneu = false;
    
    #ifdef CalNeuOscSUN_DEBUG
        printf("Configuring boxes.\n");
    #endif
    
    if(caja[0]->neutype == 0 or caja[0]->neutype == 1) {
        // one and only one neutrino type
        caja[0]->param->neutype = caja[0]->neutype;
        neu_and_aneu = false;
    } else {
        // neutrinos and antineutrinos
        neu_and_aneu = true;
    }

    // initial box setup
    if(neu_and_aneu){
        caja[0]->neutype = 0;
        caja[0]->param->neutype = 0;
        
        caja[1]->neutype = 1;
        caja[1]->param->neutype = 1;
    }

    InitSUBox(caja[0]);
    InitSUBox(caja[1]);
    
    #ifdef CalNeuOscSUN_DEBUG
        printf("Begin convert VectorState to SUState.\n");
    #endif
    
    // setup initial state    
    SUState* state_in = ConvertVectorStateToSUState(state_vector_in[0],caja[0]);

    SUState* state_in_aux = NULL;
    if (neu_and_aneu){
        state_in_aux = ConvertVectorStateToSUState(state_vector_in[1],caja[1]);
    }
    
    #ifdef CalNeuOscSUN_DEBUG
        printf("End convert VectorState to SUState.\n");
    #endif
    
    #ifdef CalNeuOscSUN_DEBUG
        printf("Convert SUState to double.\n");
    #endif

    // convert to double
    double* state_dbl = (double*) state_in;
    double* state_dbl_aux = (double*) state_in_aux;
    
    // setting up GSL ODE solver
    gsl_odeiv2_step *s = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45, caja[0]->numeqn);
    gsl_odeiv2_control *c = gsl_odeiv2_control_y_new(abs_error,rel_error);
    gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc(caja[0]->numeqn);
    
    gsl_odeiv2_step *s_aux = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45, caja[1]->numeqn);
    gsl_odeiv2_control *c_aux = gsl_odeiv2_control_y_new(abs_error,rel_error);
    gsl_odeiv2_evolve *e_aux = gsl_odeiv2_evolve_alloc(caja[1]->numeqn);
    
    // ODE system
    gsl_odeiv2_system sys = {&RHS_SUN, NULL, caja[0]->numeqn, caja[0]};
    gsl_odeiv2_system sys_aux = {&RHS_SUN, NULL, caja[1]->numeqn, caja[1]};
    
    // defining ODE extra variables
    double km = caja[0]->param->km;
    
    double x = 0;                       // ODE independent variable
    double x_aux = 0;                   // ODE2 independent variable
    double x_ini = caja[0]->track->xini;  // initial position
    double x_end = caja[0]->track->xend;  // final position
    // step sizes
    double h        = MIN(10.0*km,x_end/10.0);
    double h_min    = 1.0e-5*km;
    double h_max    = MIN(100.0*km,x_end/5.0);
    
    #ifdef USE_GSL_ODE_DRIVER
        // primary driver
        gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rkf45,abs_error,rel_error,0.0);
        gsl_odeiv2_step_set_driver(s,d);
        //gsl_odeiv2_driver_set_hmin(d,h_min);
        gsl_odeiv2_driver_set_hmax(d,h_max);
        
        // auxiliary driver
        gsl_odeiv2_driver *d_aux = gsl_odeiv2_driver_alloc_y_new(&sys_aux,gsl_odeiv2_step_rkf45,abs_error,rel_error,0.0);
        gsl_odeiv2_step_set_driver(s_aux,d_aux);
        //gsl_odeiv2_driver_set_hmin(d_aux,h_min);
        gsl_odeiv2_driver_set_hmax(d_aux,h_max);
    #endif  
    
    #ifdef CalNeuOscSUN_DEBUG
        printf("GSL paramers :\n");
        printf("x_ini : %g [km]\n", x_ini/km);
        printf("x_end : %g [km]\n", x_end/km);
        printf("h : %g [km]\n", h/km);
        printf("h_min : %g [km]\n", h_min/km);      
    #endif
    
    // initial position
    x = x_ini;
    x_aux = x_ini;
    
    #ifdef CalNeuOscSUN_DEBUG
        int count = 0;
        int count_step = 10;
    #endif
    
    #ifdef USE_FIX_STEP
        h = 1.0*km;
    #endif

    #ifdef CalNeuOscSUN_DEBUG
        printf("Start calculation.\n");
    #endif
    
    while (x < x_end){
        #ifdef USE_GSL_ODE_DRIVER
            double x_inter = x + 50.0*km;
            gsl_status = gsl_odeiv2_driver_apply(d,&x,x_inter,state_dbl);
            if( neu_and_aneu ){
                gsl_status = gsl_odeiv2_driver_apply(d_aux,&x_aux,x_inter,state_dbl_aux);
            }
        #else
            gsl_status = gsl_odeiv2_evolve_apply(e,c,s,&sys,&x,x_end,&h,state_dbl);
            if( neu_and_aneu ){
                gsl_status = gsl_odeiv2_evolve_apply(e_aux,c_aux,s_aux,&sys_aux,&x_aux,x_end,&h,state_dbl_aux);
            }
        #endif


        #ifdef ManualTauReinjection
            if(caja[0]->tau_regeneration and not neu_and_aneu){
                TauLeptonReinj_SUN(state_dbl,caja[0]);
            } else if (caja[0]->tau_regeneration)  {
                TauLeptonReinj_SUN(state_dbl,caja[0],state_dbl_aux,caja[1]);
            }
        #endif
        
        #ifdef CalNeuOscSUN_DEBUG
            if(count%count_step == 0){
                printf("x_current : %g %g [km]\n", x/km,x_aux/km);
                //printf("step : %g [km]\n", h/km);
            }
        #endif

        if( gsl_status != GSL_SUCCESS ){
            fprintf(stderr,"CalNeuOscGSL: Error in GSL ODE solver.\n");
            break;
        }
        
        if(h < h_min){
            h = h_min;
        }
        
        if(h > h_max){
            h = h_max;
        }
        
        #ifdef CalNeuOscSUN_DEBUG
            count++;
        #endif
    }
    
    #ifdef CalNeuOscSUN_DEBUG
        printf("End calculation.\n");
    #endif
    
    if (e) { gsl_odeiv2_evolve_free(e);      e = NULL; }
    if (c) { gsl_odeiv2_control_free(c);     c = NULL; }
    if (s) { gsl_odeiv2_step_free(s);        s = NULL; }
    #ifdef USE_GSL_ODE_DRIVER
    if (d) { gsl_odeiv2_driver_free(d);        d = NULL; }
    if (d_aux) { gsl_odeiv2_driver_free(d_aux);        d_aux = NULL; }    
    #endif

    // convert back to SUState
    SUState* state_out = (SUState*) state_dbl;
    SUState* state_out_aux = (SUState*) state_dbl_aux;
    
    state_vector_in[0] = ConvertSUStateToVectorState(state_out,caja[0]);
    if(neu_and_aneu){
        state_vector_in[1] = ConvertSUStateToVectorState(state_out_aux,caja[1]);
    }

    return GSL_SUCCESS;
}

int CalNeuOscSUN(VectorFluxState* state_vector_in[2], OscSetup* setup,double abs_error,double rel_error){
    int gsl_status = GSL_SUCCESS;
    bool neu_and_aneu = false;
    
    #ifdef CalNeuOscSUN_DEBUG
        printf("Configuring boxes.\n");
    #endif
    
    SUBox * caja = new SUBox;
    
    ConvertOscSetupToSUBox(setup,caja);

    if(caja->neutype == 0 or caja->neutype == 1) {
        // one and only one neutrino type
        caja->param->neutype = caja->neutype;
        neu_and_aneu = false;
    } else {
        // neutrinos and antineutrinos
        neu_and_aneu = true;
    }

    // initial box setup
    if(neu_and_aneu){
        caja->neutype = 0;
        caja->param->neutype = 0;
    }

    InitSUBox(caja);
     
   // std::cout << "lala" << neu_and_aneu << " " << caja->numeqn <<std::endl;
    // clone box object
    SUBox* caja_aux = new SUBox(*caja);
    if(neu_and_aneu){
        caja_aux->neutype = 1;
        caja_aux->param->neutype = 1;
    }

    InitSUBox(caja_aux);
    
    #ifdef CalNeuOscSUN_DEBUG
        printf("Begin convert VectorState to SUState.\n");
    #endif
    
    // setup initial state    
    SUState* state_in = ConvertVectorStateToSUState(state_vector_in[0],caja);

    SUState* state_in_aux = NULL;
    if (neu_and_aneu){
        state_in_aux = ConvertVectorStateToSUState(state_vector_in[1],caja_aux);
    }
    
    #ifdef CalNeuOscSUN_DEBUG
        printf("End convert VectorState to SUState.\n");
    #endif

    
    #ifdef CalNeuOscSUN_DEBUG
        printf("Convert SUState to double.\n");
    #endif

    // convert to double
    double* state_dbl = (double*) state_in;
    double* state_dbl_aux = (double*) state_in_aux;
    
    // setting up GSL ODE solver
    gsl_odeiv2_step *s = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45, caja->numeqn);
    gsl_odeiv2_control *c = gsl_odeiv2_control_y_new(abs_error,rel_error);
    gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc(caja->numeqn);
    
    gsl_odeiv2_step *s_aux = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45, caja_aux->numeqn);
    gsl_odeiv2_control *c_aux = gsl_odeiv2_control_y_new(abs_error,rel_error);
    gsl_odeiv2_evolve *e_aux = gsl_odeiv2_evolve_alloc(caja_aux->numeqn);
    
    // ODE system
    gsl_odeiv2_system sys = {&RHS_SUN, NULL, caja->numeqn, caja};
    gsl_odeiv2_system sys_aux = {&RHS_SUN, NULL, caja_aux->numeqn, caja_aux};
    
    // defining ODE extra variables
    double km = setup->param->km;
    
    double x = 0;                       // ODE independent variable
    double x_aux = 0;                   // ODE2 independent variable
    double x_ini = setup->track->xini;  // initial position
    double x_end = setup->track->xend;  // final position
    // step sizes
    double h        = MIN(10.0*km,x_end/10.0);
    double h_min    = 1.0e-5*km;
    double h_max    = MIN(100.0*km,x_end/5.0);
    
    #ifdef USE_GSL_ODE_DRIVER
        // primary driver
        gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rkf45,abs_error,rel_error,0.0);
        gsl_odeiv2_step_set_driver(s,d);
        //gsl_odeiv2_driver_set_hmin(d,h_min);
        gsl_odeiv2_driver_set_hmax(d,h_max);
        
        // auxiliary driver
        gsl_odeiv2_driver *d_aux = gsl_odeiv2_driver_alloc_y_new(&sys_aux,gsl_odeiv2_step_rkf45,abs_error,rel_error,0.0);
        gsl_odeiv2_step_set_driver(s_aux,d_aux);
        //gsl_odeiv2_driver_set_hmin(d_aux,h_min);
        gsl_odeiv2_driver_set_hmax(d_aux,h_max);
    #endif  
    
    #ifdef CalNeuOscSUN_DEBUG
        printf("GSL paramers :\n");
        printf("x_ini : %g [km]\n", x_ini/km);
        printf("x_end : %g [km]\n", x_end/km);
        printf("h : %g [km]\n", h/km);
        printf("h_min : %g [km]\n", h_min/km);      
    #endif
    
    // initial position
    x = x_ini;
    x_aux = x_ini;
    
    #ifdef CalNeuOscSUN_DEBUG
        int count = 0;
        int count_step = 10;
    #endif
    
    #ifdef USE_FIX_STEP
        h = 1.0*km;
    #endif

    #ifdef CalNeuOscSUN_DEBUG
        printf("Start calculation.\n");
    #endif
    
    while (x < x_end){
        #ifdef USE_GSL_ODE_DRIVER
            double x_inter = x + 50.0*km;
            gsl_status = gsl_odeiv2_driver_apply(d,&x,x_inter,state_dbl);
            if( neu_and_aneu ){
                gsl_status = gsl_odeiv2_driver_apply(d_aux,&x_aux,x_inter,state_dbl_aux);
            }
        #else
            gsl_status = gsl_odeiv2_evolve_apply(e,c,s,&sys,&x,x_end,&h,state_dbl);
            if( neu_and_aneu ){
                gsl_status = gsl_odeiv2_evolve_apply(e_aux,c_aux,s_aux,&sys_aux,&x_aux,x_end,&h,state_dbl_aux);
            }
        #endif

        #ifdef ManualTauReinjection
            if(setup->tau_regeneration and not neu_and_aneu){
                TauLeptonReinj_SUN(state_dbl,caja);
            } else if (setup->tau_regeneration) {
                TauLeptonReinj_SUN(state_dbl,caja,state_dbl_aux,caja_aux);
            }
        #endif
        
        #ifdef CalNeuOscSUN_DEBUG
            if(count%count_step == 0){
                printf("x_current : %g %g [km]\n", x/km,x_aux/km);
                //printf("step : %g [km]\n", h/km);
            }
        #endif

        if( gsl_status != GSL_SUCCESS ){
            fprintf(stderr,"CalNeuOscGSL: Error in GSL ODE solver.\n");
            break;
        }
        
        if(h < h_min){
            h = h_min;
        }
        
        if(h > h_max){
            h = h_max;
        }
        
        #ifdef CalNeuOscSUN_DEBUG
            count++;
        #endif
    }
    
    #ifdef CalNeuOscSUN_DEBUG
        printf("End calculation.\n");
    #endif
    
    if (e) { gsl_odeiv2_evolve_free(e);      e = NULL; }
    if (c) { gsl_odeiv2_control_free(c);     c = NULL; }
    if (s) { gsl_odeiv2_step_free(s);        s = NULL; }
    #ifdef USE_GSL_ODE_DRIVER
    if (d) { gsl_odeiv2_driver_free(d);        d = NULL; }
    if (d_aux) { gsl_odeiv2_driver_free(d_aux);        d_aux = NULL; }    
    #endif

    // convert back to SUState
    SUState* state_out = (SUState*) state_dbl;
    SUState* state_out_aux = (SUState*) state_dbl_aux;

    /*
    std::cout << state_out->neu_array[0] << std::endl;
    std::cout << "proyectors : " << std::endl;
    std::cout << caja->evol_flavor_proj[0] << std::endl;
    std::cout << caja->evol_flavor_proj[1] << std::endl;
    std::cout << caja->evol_flavor_proj[2] << std::endl;
*/
    /*
    // convert to vector state
    VectorFluxState state_vector_out[2];
    
    state_vector_out[0] = ConvertSUStateToVectorState(state_out,caja);
    if(neu_and_aneu){
        state_vector_out[1] = ConvertSUStateToVectorState(state_out_aux,caja_aux);
    }
    */
    
    state_vector_in[0] = ConvertSUStateToVectorState(state_out,caja);
    if(neu_and_aneu){
        state_vector_in[1] = ConvertSUStateToVectorState(state_out_aux,caja_aux);
    }

    //return state_vector_out;
    return GSL_SUCCESS;
};

