#include "neurho.h"

// Macros
#define SQR(x)      ((x)*(x))                        // x^2
#define SQR_ABS(x)  (SQR(creal(x)) + SQR(cimag(x)))  // |x|^2
#define POW10(x)    (exp(M_LN10*(x)))                // 10^x
#define MIN(X,Y)    ( ((X) < (Y)) ? (X) : (Y) )
#define MAX(X,Y)    ( ((X) > (Y)) ? (X) : (Y) )
#define SIGN(a,b)   ( (b) > 0.0 ? (fabs(a)) : (-fabs(a)) )
#define KRONECKER(i,j)  ( (i)==(j) ? 1 : 0 )

gsl_matrix_complex* Projector(int flavor,int numneu){
    gsl_matrix_complex* proj = gsl_matrix_complex_calloc(numneu,numneu);
    gsl_matrix_complex_set(proj,flavor,flavor,gsl_complex_rect(1.0,0.0));
    return proj;
}

gsl_matrix_complex* FlavorProjector(int flavor,int numneu){
    gsl_matrix_complex* proj = gsl_matrix_complex_calloc(numneu,numneu);
    gsl_matrix_complex_set(proj,flavor,flavor,gsl_complex_rect(1.0,0.0));
    return proj;
}

gsl_matrix_complex* FlavorProjector(int flavor,OscSetup* setup){
    int numneu = setup->param->numneu;
    gsl_matrix_complex* proj = gsl_matrix_complex_calloc(numneu,numneu);
    
    gsl_matrix_complex_set(proj,flavor,flavor,gsl_complex_rect(1.0,0.0));
    
    if(setup->basis == 1){
        // if on mass basis
        gsl_matrix_complex *U1 = gsl_matrix_complex_alloc(numneu,numneu);
        gsl_matrix_complex *U2 = gsl_matrix_complex_alloc(numneu,numneu);
        // temporal matrices to make the multiplication
        gsl_matrix_complex *T1 = gsl_matrix_complex_alloc(numneu,numneu);
        
        gsl_matrix_complex_memcpy(U1,setup->U);
        gsl_matrix_complex_memcpy(U2,setup->U);
        
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,
                       gsl_complex_rect(1.0,0.0),proj,
                       U1,gsl_complex_rect(0.0,0.0),T1);
        gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,
                       gsl_complex_rect(1.0,0.0),U2,
                       T1,gsl_complex_rect(0.0,0.0),proj);
        
        gsl_matrix_complex_free(U1);
        gsl_matrix_complex_free(U2);
        gsl_matrix_complex_free(T1);
        
    } else if(setup->basis == 2){
        // if on interaction basis use mass basis and update at every step
        gsl_matrix_complex *U1 = gsl_matrix_complex_alloc(numneu,numneu);
        gsl_matrix_complex *U2 = gsl_matrix_complex_alloc(numneu,numneu);
        // temporal matrices to make the multiplication
        gsl_matrix_complex *T1 = gsl_matrix_complex_alloc(numneu,numneu);
        
        gsl_matrix_complex_memcpy(U1,setup->U);
        gsl_matrix_complex_memcpy(U2,setup->U);
        
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,
                       gsl_complex_rect(1.0,0.0),proj,
                       U1,gsl_complex_rect(0.0,0.0),T1);
        gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,
                       gsl_complex_rect(1.0,0.0),U2,
                       T1,gsl_complex_rect(0.0,0.0),proj);
        
        gsl_matrix_complex_free(U1);
        gsl_matrix_complex_free(U2);
        gsl_matrix_complex_free(T1);
        
    }

    return proj;
}

vector<gsl_matrix_complex*> InitProjectorVector(OscSetup* setup){
    int numneu = setup->param->numneu;
    vector<gsl_matrix_complex*> proyectors;
    for(int i = 0; i < numneu; i++){
        proyectors.push_back(FlavorProjector(i,setup));
    }
    return proyectors;
}

gsl_vector* TraceFlavorProject(int flavor,VectorFluxState* state_vector){
    FluxState* state_zero = state_vector->at(0);
    int numneu = state_zero->F_nu->size1;
    
    gsl_vector* result = gsl_vector_calloc(state_vector->size());
    
    for(unsigned int e = 0; e < state_vector->size(); e++){
        FluxState* state = state_vector->at(e);
        gsl_matrix_complex* product = gsl_matrix_complex_calloc(numneu,numneu);
        gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1.0,0.0),
                       FlavorProjector(flavor,numneu),
                       state->F_nu,
                       gsl_complex_rect(0.0,0.0),product);
        
        double trace = 0.0;
        for(int i = 0; i<numneu; i++){
            trace = trace + gsl_matrix_complex_get(product,i,i).dat[0];
        }
        gsl_vector_set(result,e,trace);
    }
    
    return result;
}

InterOperator* InitInter(OscSetup* setup){
    int e_size = setup->E_range.size();
    PhysConst* param = setup->param;
    int neutype = setup->neutype;
    // initialize vectors
    InterOperator* inter = new InterOperator;

    // load cross sections
    NeutrinoCrossSections ncs(setup->E_range[0],setup->E_range[e_size-1],e_size-1);
    
    // initializing cross section pointers
    for(int flavor = 0; flavor < setup->param->numneu; flavor++){
        // differential cross sections matrices
        gsl_matrix_complex* dsignudE_NC = gsl_matrix_complex_calloc(e_size,e_size);
        gsl_matrix_complex* dsignudE_CC = gsl_matrix_complex_calloc(e_size,e_size);
        
        gsl_matrix_complex* dNledE_NC = gsl_matrix_complex_calloc(e_size,e_size);
        gsl_matrix_complex* dNledE_CC = gsl_matrix_complex_calloc(e_size,e_size);
        
        // total cross sections arrays
        gsl_vector_complex* sigma_CC = gsl_vector_complex_calloc(e_size);
        gsl_vector_complex* sigma_NC = gsl_vector_complex_calloc(e_size);

        for(int e1 = 0; e1 < e_size; e1++){
            // differential cross sections
            for(int e2 = 0; e2 < e1; e2++){
                // neutrino
                gsl_matrix_complex_set(dsignudE_NC,e1,e2,
                                       gsl_complex_rect(pow(param->cm,2)*pow(param->GeV,-1)*
                                                        ncs.dsde_NC(e1,e2,flavor,neutype),
                                                        0.0)
                                       );

                gsl_matrix_complex_set(dsignudE_CC,e1,e2,
                                       gsl_complex_rect(pow(param->cm,2)*pow(param->GeV,-1)*
                                                        ncs.dsde_CC(e1,e2,flavor,neutype),
                                                        0.0)
                                       );
            }
            // total cross sections
            gsl_vector_complex_set(sigma_CC,e1,
                                   gsl_complex_rect(pow(param->cm,2)*
                                                    ncs.sigma_CC(e1,flavor,neutype),
                                                    0.0)
                                   );
            gsl_vector_complex_set(sigma_NC,e1,
                                   gsl_complex_rect(pow(param->cm,2)*
                                                    ncs.sigma_NC(e1,flavor,neutype),
                                                    0.0)
                                   );
        }
        
        
        #ifdef InitInter_Unitarity_FIX
            // cc cross rescaling
            gsl_complex sig_cc_min = gsl_vector_complex_get(sigma_CC,0);
            
            for(int e1 = 0; e1 < e_size; e1++){
                gsl_complex sig_cc_int = gsl_complex_rect(0.0,0.0);
                gsl_complex sig_cc = gsl_vector_complex_get(sigma_CC,e1);
                
                for(int e2 = 0; e2 < e1; e2++){
                    gsl_complex sig_cc_int_term = gsl_complex_mul_real(
                                                        gsl_matrix_complex_get(dsignudE_CC,e1,e2),
                                                        (setup->E_range[e2]-setup->E_range[e2-1])
                                                        );
                    
                    sig_cc_int = gsl_complex_add(sig_cc_int,sig_cc_int_term);
                }
                
                if(e1 != 0){
                    gsl_complex int_sig_rescale = gsl_complex_div(gsl_complex_sub(sig_cc,sig_cc_min),sig_cc_int);
                    
                    for(int e2 = 0; e2 < e1; e2++){
                        gsl_complex dsde_cc = gsl_matrix_complex_get(dsignudE_CC,e1,e2);
                        gsl_matrix_complex_set(dsignudE_CC,e1,e2,gsl_complex_mul(int_sig_rescale,dsde_cc));
                    }
                }
                
                #ifdef InitInter_Unitarity_DEBUG
                    // recheck integral
                    gsl_complex sig_cc_int_check = gsl_complex_rect(0.0,0.0);
                    for(int e2 = 0; e2 < e1; e2++){
                        gsl_complex sig_cc_int_check_term = gsl_complex_mul_real(
                                                            gsl_matrix_complex_get(dsignudE_CC,e1,e2),
                                                            (setup->E_range[e2]-setup->E_range[e2-1])
                                                            );
                        
                        sig_cc_int_check = gsl_complex_add(sig_cc_int_check_term,sig_cc_int_check);
                    }
                if(e1 != 0){
                    cout << "== dsde_CC compare ==" << endl;
                    double sig = gsl_vector_complex_get(sigma_CC,e1).dat[0];
                    double sig_int = sig_cc_int.dat[0];
                    double sig_res = sig_cc_int_check.dat[0];
                    double sig_min = sig_cc_min.dat[0];
                    cout << "E : " << setup->E_range[e1]*1.0e-9 << " sig : " << sig<<  " sig_int " << sig_int <<
                    " sig_res " << sig_res << " r " << sig_int/sig << " r_res " << sig_res/(sig-sig_min)<< endl;
                }
                #endif
            }
            
            // nc cross rescaling
            gsl_complex sig_nc_min = gsl_vector_complex_get(sigma_NC,0);
            
            for(int e1 = 0; e1 < e_size; e1++){
                gsl_complex sig_nc_int = gsl_complex_rect(0.0,0.0);
                gsl_complex sig_nc = gsl_vector_complex_get(sigma_NC,e1);
                
                for(int e2 = 0; e2 < e1; e2++){
                    gsl_complex sig_nc_int_term = gsl_complex_mul_real(
                                                        gsl_matrix_complex_get(dsignudE_NC,e1,e2),
                                                        (setup->E_range[e2]-setup->E_range[e2-1])
                                                        );
                    
                    sig_nc_int = gsl_complex_add(sig_nc_int,sig_nc_int_term);
                }
                
                if(e1 != 0){
                    gsl_complex int_sig_rescale = gsl_complex_div(gsl_complex_sub(sig_nc,sig_nc_min),sig_nc_int);
                    
                    for(int e2 = 0; e2 < e1; e2++){
                        gsl_complex dsde_nc = gsl_matrix_complex_get(dsignudE_NC,e1,e2);
                        gsl_matrix_complex_set(dsignudE_NC,e1,e2,gsl_complex_mul(int_sig_rescale,dsde_nc));
                    }
                }
                
                #ifdef InitInter_Unitarity_DEBUG
                    // recheck integral
                    gsl_complex sig_nc_int_check = gsl_complex_rect(0.0,0.0);
                    for(int e2 = 0; e2 < e1; e2++){
                        gsl_complex sig_nc_int_check_term = gsl_complex_mul_real(
                                                            gsl_matrix_complex_get(dsignudE_NC,e1,e2),
                                                            (setup->E_range[e2]-setup->E_range[e2-1])
                                                            );
                        
                        sig_nc_int_check = gsl_complex_add(sig_nc_int_check_term,sig_nc_int_check);
                    }
                if(e1 != 0){
                    cout << "== dsde_NC compare ==" << endl;
                    double sig = gsl_vector_complex_get(sigma_NC,e1).dat[0];
                    double sig_int = sig_nc_int.dat[0];
                    double sig_res = sig_nc_int_check.dat[0];
                    double sig_min = sig_nc_min.dat[0];
                    cout << "E : " << setup->E_range[e1]*1.0e-9 << " sig : " << sig<<  " sig_int " << sig_int <<
                    " sig_res " << sig_res << " r " << sig_int/sig << " r_res " << sig_res/(sig-sig_min)<< endl;
                }
                #endif
            }
        #endif

        // push to structure for differential cross sections
        inter->dsignudE_NC.push_back(dsignudE_NC);
        inter->dsignudE_CC.push_back(dsignudE_CC);
        
        inter->dNledE_NC.push_back(dNledE_NC);
        inter->dNledE_CC.push_back(dNledE_CC);
        // push to structure for total cross sections
        inter->sigma_CC.push_back(sigma_CC);
        inter->sigma_NC.push_back(sigma_NC);
    }
    
    // initializing interaction length pointers
    for(int flavor = 0; flavor < setup->param->numneu; flavor++){
        gsl_vector_complex* L_int = gsl_vector_complex_calloc(e_size);
        gsl_vector_complex* L_NC = gsl_vector_complex_calloc(e_size);
        gsl_vector_complex* L_CC = gsl_vector_complex_calloc(e_size);
        
        inter->Lint.push_back(L_int);
        inter->LNC.push_back(L_NC);
        inter->LCC.push_back(L_CC);
    }
    
    // initializing tau decay length pointer
    inter->LTauDecay = gsl_vector_complex_calloc(e_size);
    
    for(int e = 0; e < e_size; e++){
        gsl_vector_complex_set(inter->LTauDecay,e,
                               gsl_complex_rect(setup->param->tau_lifetime*
                                                setup->E_range[e]/
                                                setup->param->tau_mass
                                                ,0.0));
    }
    
    // load tau decay spectra
    TauDecaySpectra tdc(setup->E_range[0],setup->E_range[e_size-1],e_size-1);
    
    // initializing tau decay spectra
    gsl_matrix_complex* dNdENuTau_All = gsl_matrix_complex_calloc(e_size,e_size);
    gsl_matrix_complex* dNdENuTau_Lep = gsl_matrix_complex_calloc(e_size,e_size);
    
    gsl_matrix_complex* dNdELeTau_All = gsl_matrix_complex_calloc(e_size,e_size);
    gsl_matrix_complex* dNdELeTau_Lep = gsl_matrix_complex_calloc(e_size,e_size);
    
    for(int e1 = 0; e1 < e_size; e1++){
        for(int e2 = 0; e2 < e1; e2++){
            gsl_matrix_complex_set(dNdENuTau_All,e1,e2,
                                       gsl_complex_rect(pow(param->GeV,-1)*
                                                        tdc.dNdEnu_All(e1,e2),
                                                        0.0)
                                       );
            
            gsl_matrix_complex_set(dNdENuTau_Lep,e1,e2,
                                       gsl_complex_rect(pow(param->GeV,-1)*
                                                        tdc.dNdEnu_Lep(e1,e2),
                                                        0.0)
                                       );
            
            gsl_matrix_complex_set(dNdELeTau_All,e1,e2,
                                       gsl_complex_rect(pow(param->GeV,-1)*
                                                        tdc.dNdEle_All(e1,e2),
                                                        0.0)
                                       );
            
            gsl_matrix_complex_set(dNdELeTau_Lep,e1,e2,
                                       gsl_complex_rect(pow(param->GeV,-1)*
                                                        tdc.dNdEle_Lep(e1,e2),
                                                        0.0)
                                       );
        }
    }
    
    inter->dNdENuTau_All = dNdENuTau_All;
    inter->dNdENuTau_Lep = dNdENuTau_Lep;
    inter->dNdELeTau_All = dNdELeTau_All;
    inter->dNdELeTau_Lep = dNdELeTau_Lep;
    
    #ifdef InitInter_DEBUG
        cout << "== TauDecaySpectra - all ==" << endl;
        for(int e1 = 0; e1 < e_size; e1++){]
            for(int e2 = 0; e2 < e_size; e2++){
                cout << gsl_matrix_complex_get(dNdENuTau_All,e1,e2).dat[0] <<
                "+i" << gsl_matrix_complex_get(dNdENuTau_All,e1,e2).dat[1] << " ";
            }
            cout << endl;
        }
    #endif
    
    return inter;
}

int UpdateInter(InterOperator* inter,OscSetup* setup){
    int e_size = setup->E_range.size();
    //double ye = setup->body->ye(setup->track);
    double density = setup->body->density(setup->track);
    //double dens = density*setup->param->gr*pow(setup->param->cm,-3);

    PhysConst* param = setup->param;
    double num_nuc = param->Na*pow(param->cm,-3)*density;//*(1.0-ye)/2.0;
    if (num_nuc == 0){
        num_nuc = param->Na*pow(param->cm,-3)*1.0e-10;
    }

    #ifdef UpdateInter_DEBUG
        cout << "== INTER at x = " << (setup->track->x)/(setup->param->km) << " [km] =="<< endl;
        cout << "NumNuc : "<< num_nuc << endl;
        cout << "==" << endl;
    #endif

    // calculating charge current/neutral current lengths
    
    for(int flavor = 0; flavor < MIN(param->numneu,3); flavor++){
        #ifdef UpdateInter_DEBUG
            cout << "============" << flavor << "============" << endl;
        #endif 
        // load interaction length vectors pointers
        gsl_vector_complex* Lint = inter->Lint[flavor];
        gsl_vector_complex* LNC = inter->LNC[flavor];
        gsl_vector_complex* LCC = inter->LCC[flavor];
        // load cross section pointers
        gsl_vector_complex* sigNC = inter->sigma_NC[flavor];
        gsl_vector_complex* sigCC = inter->sigma_CC[flavor];
        
        for(int e = 0; e < e_size ; e++){
            
            gsl_complex NC = gsl_complex_inverse(
                                        gsl_complex_mul_real(
                                            gsl_vector_complex_get(sigNC,e),
                                            num_nuc)
                                        );
            gsl_complex CC = gsl_complex_inverse(
                                        gsl_complex_mul_real(
                                            gsl_vector_complex_get(sigCC,e),
                                            num_nuc)
                                        );
            
            #ifdef UpdateInter_DEBUG
            cout << "== CC NC Terms x = " << (setup->track->x)/(setup->param->km)
            << " [km]  E = " << setup->E_range[e] << " [eV] =="<< endl;
            cout << "CC : " << CC.dat[0] << " NC : " << NC.dat[0] << endl;
            cout << "==" << endl;
            #endif
            
            gsl_complex NC_CC = gsl_complex_div(gsl_complex_mul(NC,CC),gsl_complex_add(NC,CC));
            
            gsl_vector_complex_set(LNC,e,NC);
            gsl_vector_complex_set(LCC,e,CC);
            gsl_vector_complex_set(Lint,e,NC_CC);
        }
    }
    
    // calculate tay decay length
    for(int e = 0; e < e_size; e++){
        gsl_vector_complex_set(inter->LTauDecay,e,
                                gsl_complex_rect(setup->param->tau_lifetime*
                                                setup->E_range[e]/
                                                setup->param->tau_mass
                                                ,0.0));
    }
    
    #ifdef UpdateInter_DEBUG
        cout << "== TauLength Terms x = " << (setup->track->x)/(setup->param->km) << " [km] ==" <<endl;
        for(int e = 0; e < e_size; e++){
            cout << gsl_vector_complex_get(inter->LTauDecay,e).dat[0] << " ";
        }
        cout << endl;
        cout << "==" << endl;
    #endif
    
    return GSL_SUCCESS;
}

VectorFluxState* InitSimpleState(int ini_flavor,OscSetup *setup){
    VectorFluxState* state_vector  = new VectorFluxState;
    
    int e_size = setup->E_range.size();
    // create initial state : neutrino density matrix
    int numneu = setup->param->numneu;
    for(int e=0; e<e_size; e++){
        // initialize neutrino flux density matrix
        gsl_matrix_complex* F_nu = gsl_matrix_complex_calloc(numneu,numneu);
        gsl_matrix_complex_set(F_nu,
                                   ini_flavor,ini_flavor,
                                   gsl_complex_rect(1,0));
        // initialize lepton matrix
        gsl_matrix_complex* F_le = gsl_matrix_complex_calloc(numneu,1);
        // putting into container
        FluxState* state = new FluxState;
        state->F_nu = F_nu;
        state->F_le = F_le;
        // putting into state vector
        state_vector->push_back(state);
    }
    
    return state_vector;
}

VectorFluxState* InitPowerLawFluxState(int ini_flavor,double ene_power, OscSetup *setup){
    VectorFluxState* state_vector  = new VectorFluxState;
    
    int e_size = setup->E_range.size();
    // create initial state : neutrino density matrix
    int numneu = setup->param->numneu;
    for(int e=0; e<e_size; e++){
        // initialize neutrino flux density matrix
        gsl_matrix_complex* F_nu = gsl_matrix_complex_calloc(numneu,numneu);
        if (ene_power == 0)
          gsl_matrix_complex_set(F_nu,
                                     ini_flavor,ini_flavor,
                                     gsl_complex_rect(1.0,0.0));
        else
          gsl_matrix_complex_set(F_nu,
                                     ini_flavor,ini_flavor,
                                     gsl_complex_rect(pow(setup->E_range[e],ene_power),0));
        // initialize lepton matrix
        gsl_matrix_complex* F_le = gsl_matrix_complex_calloc(numneu,1);
        // putting into container
        FluxState* state = new FluxState;
        state->F_nu = F_nu;
        state->F_le = F_le;
        // putting into state vector
        state_vector->push_back(state);
    }
    return state_vector;
}

VectorFluxState* InitBodyFluxState(OscSetup* setup){
    VectorFluxState* state_vector  = new VectorFluxState;
    
    if (setup->body->name == "SunASnu"){
        // load some data
        Table nu_KL_flux_tbl = quickread("/Users/carguelles/Workspace/atm_solar_neu/data/flux/K_L.txt");
        Table nu_Kpm_flux_tbl = quickread("/Users/carguelles/Workspace/atm_solar_neu/data/flux/K_plus_minus.txt");
        Table nu_mupm_flux_tbl = quickread("/Users/carguelles/Workspace/atm_solar_neu/data/flux/mu_plus_minus.txt");
        Table nu_pionpm_flux_tbl = quickread("/Users/carguelles/Workspace/atm_solar_neu/data/flux/pi_plus_minus.txt");
        // number of points
        int e_size = setup->E_range.size();
        // create GeV Energy array
        Row E_range_GeV;
        for(int i=0; i<e_size; i++){
            E_range_GeV.push_back(setup->E_range[i]/setup->param->GeV);
        }
        // calculate new tables for the given energy range
        Table nu_KL_tbl = intertable(nu_KL_flux_tbl,E_range_GeV,0,1);
        Table nu_Kpm_tbl = intertable(nu_Kpm_flux_tbl,E_range_GeV,0,1);
        Table nu_mupm_tbl = intertable(nu_mupm_flux_tbl,E_range_GeV,0,1);
        Table nu_pionpm_tbl = intertable(nu_pionpm_flux_tbl,E_range_GeV,0,1);
        
        // calculate the combined flux
        Table nu_flux_tbl;
        for(int i = 0; i<e_size; i++){
            Row row;
            row.push_back(setup->E_range[i]);
            row.push_back(nu_KL_tbl[i][1]+nu_Kpm_tbl[i][1]+nu_mupm_tbl[i][1]+nu_pionpm_tbl[i][1]);
            
            nu_flux_tbl.push_back(row);
        }
        
        // create initial state : neutrino flux density matrix
        int numneu = setup->param->numneu;
        for(int e=0; e<e_size; e++){
            // initialize neutrino flux density matrix
            gsl_matrix_complex* F_nu = gsl_matrix_complex_calloc(numneu,numneu);
            gsl_matrix_complex_set(F_nu,
                                   setup->param->muon,setup->param->muon,
                                   gsl_complex_rect(nu_flux_tbl[e][1],0));
            // initialize lepton matrix
            gsl_matrix_complex* F_le = gsl_matrix_complex_calloc(numneu,1);
            // putting into container
            FluxState* state = new FluxState;
            state->F_nu = F_nu;
            state->F_le = F_le;
            // putting into state vector
            state_vector->push_back(state);
        } 
    } else if (setup->body->name == "EarthAtm"){
        int e_size = setup->E_range.size();
        int numneu = setup->param->numneu;
        for(int e=0; e<e_size; e++){
            // initialize neutrino flux density matrix
            gsl_matrix_complex* F_nu = gsl_matrix_complex_calloc(numneu,numneu);
            gsl_matrix_complex_set(F_nu,
                                   setup->param->muon,setup->param->muon,
                                   gsl_complex_rect(1.0/setup->E_range[e],0));
            // initialize lepton matrix
            gsl_matrix_complex* F_le = gsl_matrix_complex_calloc(numneu,1);
            // putting into container
            FluxState* state = new FluxState;
            state->F_nu = F_nu;
            state->F_le = F_le;
            // putting into state vector
            state_vector->push_back(state);
        }
    }
    return state_vector;
}

gsl_vector_complex* ConvertVectorStateToVector(VectorFluxState* state_vector,OscSetup* setup){
    int numneu = setup->param->numneu;
    int state_vector_size = state_vector->size();
    int complex_vector_size = state_vector_size*(SQR(numneu)+numneu);

    #ifdef ConvertVectorStateToVector_DEBUG
        cout << "numneu : " << numneu << endl;
        cout << "state vector size : " << state_vector_size << endl;
        cout << "complex vector size : " << complex_vector_size << endl;
    #endif
    
    gsl_vector_complex* vector = gsl_vector_complex_calloc(complex_vector_size);
    for(int e = 0; e < state_vector_size; e++){
        FluxState* state = state_vector->at(e);
        
        if (setup-> basis == 1){
            // go to mass matrix
            gsl_matrix_complex_change_basis_UCMU(setup->U,state->F_nu);
        } else if (setup-> basis == 2){
            // go to mass matrix
            gsl_matrix_complex_change_basis_UCMU(setup->U,state->F_nu);
            // done since mass basis = interaction basis at t = 0
        }
        
        for(int i = 0; i < SQR(numneu)+numneu; i ++){
            if(i < SQR(numneu)){
                // neutrino part
                #ifdef ConvertVectorStateToVector_DEBUG
                    cout << "vector_entry : " << e*SQR(numneu)+i << endl;
                    gsl_complex z = gsl_matrix_complex_get(state->F_nu,i/numneu,i%numneu);
                    cout << "vector_value : " << z.dat[0] << "+i" << z.dat[1] << endl;
                #endif
                gsl_vector_complex_set(vector,e*SQR(numneu)+i,
                                        gsl_matrix_complex_get(state->F_nu,i/numneu,i%numneu));
            } else {
                // lepton part
                #ifdef ConvertVectorStateToVector_DEBUG
                    cout << "vector_entry : " << SQR(numneu)*state_vector_size+e*numneu+(i-SQR(numneu)) << endl;
                    gsl_complex z = gsl_matrix_complex_get(state->F_le,i-SQR(numneu),0);
                    cout << "vector_value : " << z.dat[0] << "+i" << z.dat[1] << endl;
                #endif      
                gsl_vector_complex_set(vector,SQR(numneu)*state_vector_size+e*numneu+(i-SQR(numneu)),
                                       gsl_matrix_complex_get(state->F_le,i-SQR(numneu),0));
            }
        }
        free(state);
    }
    return vector;
}

VectorFluxState* ConvertVectorToVectorState(gsl_vector_complex* complex_vector,OscSetup* setup){
    int numneu = setup->param->numneu;
    int e_size = (complex_vector->size/(SQR(numneu)+numneu));
    
    VectorFluxState* state_vector = new VectorFluxState;
    for(int e = 0; e < e_size; e++){
        FluxState* state = new FluxState;
        // initialize neutrino flux density matrix
        state->F_nu = gsl_matrix_complex_calloc(numneu,numneu);
        for(int i = 0; i < numneu; i++){
            for(int j = 0; j < numneu; j++){
                gsl_matrix_complex_set(state->F_nu,i,j,
                                       gsl_vector_complex_get(complex_vector,e*SQR(numneu)+i*numneu+j));
            }
        }
        
        if (setup-> basis == 1){
            // go back to flavor basis from mass basis
            gsl_matrix_complex_change_basis_UMUC(setup->U,state->F_nu);
        } else if (setup-> basis == 2){
            // go to mass matrix from interaction basis
            gsl_matrix_complex* Uint = massExpH0(setup->param,setup->E_range[e],
                                                setup->track,setup->massM2);
            gsl_matrix_complex_change_basis_UCMU(Uint,state->F_nu);
            gsl_matrix_complex_free(Uint);
            // go back to flavor basis from mass basis
            gsl_matrix_complex_change_basis_UMUC(setup->U,state->F_nu);
        }

        // initialize lepton matrix
        state->F_le = gsl_matrix_complex_calloc(numneu,1);
        for(int i = 0; i < numneu; i++){
            gsl_matrix_complex_set(state->F_le,i,0,
                                   gsl_vector_complex_get(complex_vector,e_size*SQR(numneu)+e*numneu+i));
        }
        //push back
        state_vector->push_back(state);
    }
    
    return state_vector;
}

int TauLeptonReinjection(double *state_dbl, void *box){
    // Format container
    Box* caja = (Box*) box;
    // init variables
    vector<double> E_range = caja->setup->E_range;
    int numneu = caja->setup->param->numneu;
    //int numstate = SQR(numneu)+numneu;
    int e_size = E_range.size();
    
    int lep_tau = caja->setup->param->tau;
    int neu_tau = numneu*lep_tau+lep_tau;
    
    // convert to complex state
    gsl_vector_complex* state = gsl_vector_complex_calloc(caja->numeqn/2);
    state->data = (double*) state_dbl;
    
    // load the cross sections and decay spectra
    gsl_matrix_complex* dNdENuTau_All = caja->interactions->dNdENuTau_All;
    
    #ifdef TauLeptonReinjection_DEBUG
    // total tau flux integral
    gsl_complex leptau_integral = gsl_complex_rect(0.0,0.0);
    for(int e1 = 0; e1 < e_size; e1++){
        gsl_complex leptau_integral_term = gsl_complex_mul_real(
                                                        gsl_vector_complex_get(state,SQR(numneu)*e_size + numneu*e1 + lep_tau),
                                                        (E_range[e1]-E_range[e1-1])
                                                        );
        leptau_integral = gsl_complex_add(leptau_integral,leptau_integral_term);
    }
    cout << "== Lepton Integral ==" << endl;
    cout << leptau_integral.dat[0] << " + i " << leptau_integral.dat[1] << endl;
    #endif
    
    #ifdef TauLeptonReinjection_DEBUG
    // checking unitarity of dNdENuTau_All under this aproximation
    gsl_vector_complex* taudecay_integral_vector = gsl_vector_complex_calloc(e_size);
    for(int e1 = 0; e1 < e_size; e1++){
        gsl_complex taudecay_integral = gsl_complex_rect(0.0,0.0);
        for(int e2 = e1 + 1; e2 < e_size; e2++){
            gsl_complex taudecay_integral_term = gsl_complex_mul_real(
                                                            gsl_matrix_complex_get(dNdENuTau_All,e2,e1),
                                                            (E_range[e2]-E_range[e2-1])
                                                            );
            taudecay_integral = gsl_complex_add(taudecay_integral,taudecay_integral_term);
            
            //cout << E_range[e1]*1.0e-9 << " " << E_range[e2]*1.0e-9 << " " << gsl_matrix_complex_get(dNdELeTau_All,e2,e1).dat[0] << endl;
        }
        
        gsl_vector_complex_set(taudecay_integral_vector,e1,taudecay_integral);
        //cout << "== TauDecaySpectra Integral ==" << endl;
        //cout << "E : "<< E_range[e1]*1.0e-9 << " " << taudecay_integral.dat[0] << " + i " << taudecay_integral.dat[1] << endl;
    }
    // note that the first two bins of this are always zero. We neglect this correction.
    gsl_vector_complex_set(taudecay_integral_vector,0,gsl_complex_rect(1.0,0.0));
    gsl_vector_complex_set(taudecay_integral_vector,1,gsl_complex_rect(1.0,0.0));
    #endif
    
    #ifdef TauLeptonReinjection_DEBUG
        gsl_vector_complex* tau_neu_flux_vect = gsl_vector_complex_calloc(e_size);
    #endif
    
    for(int e1 = 0; e1 < e_size; e1++){
        gsl_complex reinj_tau_neu_flux = gsl_complex_rect(0.0,0.0);
        gsl_complex reinj_tau_neu_flux_lep = gsl_complex_rect(0.0,0.0);
        // calculate the integral 
        for(int e2 = e1 + 1 ; e2 < e_size; e2++){
            // all decays
            gsl_complex nutau_integral_term = gsl_complex_mul_real(
                                            gsl_complex_mul(
                                                gsl_matrix_complex_get(dNdENuTau_All,e2,e1),
                                                gsl_vector_complex_get(state,SQR(numneu)*e_size + numneu*e2 + lep_tau)
                                                ),
                                            (E_range[e2]-E_range[e2-1])
                                            );                    
            reinj_tau_neu_flux = gsl_complex_add(nutau_integral_term,reinj_tau_neu_flux);
        }
                
        // rescale to integrated differential cross section
        
        // reinj_tau_neu_flux = gsl_complex_div(reinj_tau_neu_flux,gsl_complex_rect(1.0,0));
        
        #ifdef TauLeptonReinjection_DEBUG
        gsl_vector_complex_set(tau_neu_flux_vect,e1,reinj_tau_neu_flux);
        #endif
    
        
        // add to current tau neutrino flux
        gsl_complex tau_neu_flux = gsl_vector_complex_get(state,SQR(numneu)*e1+neu_tau);
        tau_neu_flux = gsl_complex_add(tau_neu_flux,reinj_tau_neu_flux);
        gsl_vector_complex_set(state,SQR(numneu)*e1+neu_tau,tau_neu_flux);
    }
    
    #ifdef TauLeptonReinjection_DEBUG
    gsl_complex nutau_integral = gsl_complex_rect(0.0,0.0);
    for(int e1 = 0; e1 < e_size; e1++){
        gsl_complex nutau_integral_term = gsl_complex_mul_real(
                                                        gsl_vector_complex_get(tau_neu_flux_vect,e1),
                                                        (E_range[e1]-E_range[e1-1])
                                                        );
        nutau_integral = gsl_complex_add(nutau_integral,nutau_integral_term);
    }   
    
    cout << "== Neutrino Integral ==" << endl;
    cout << nutau_integral.dat[0] << " + i " << nutau_integral.dat[1] << endl;
    cout << "== Integral Ratio nu/lep ==" << endl;
    cout << nutau_integral.dat[0]/leptau_integral.dat[0] << endl; 
    
    gsl_vector_complex_free(tau_neu_flux_vect);
    #endif
    
    
    
    // remove all the tau lepton flux
    for(int e1 = 0; e1 < e_size; e1++){
        gsl_vector_complex_set(state,SQR(numneu)*e_size + numneu*e1 + lep_tau,gsl_complex_rect(0.0,0.0));
    }
    
    // write back to double representation state
    for (int i = 0 ; i< caja->numeqn; i++){
        state_dbl[i] = state->data[i];     
    }
    
    gsl_vector_complex_free(state);
    
    return GSL_SUCCESS;
}

int TauLeptonReinjection(double *state_dbl, void *box, double *state_dbl_aux, void *box_aux){
    // Format container
    Box* caja = (Box*) box;
    Box* caja_aux = (Box*) box_aux;
    // init variables
    vector<double> E_range = caja->setup->E_range;
    int numneu = caja->setup->param->numneu;
    //int numstate = SQR(numneu)+numneu;
    int e_size = E_range.size();
    
    // lepton indexes
    int lep_tau = caja->setup->param->tau;
    int lep_muon = caja->setup->param->muon;
    int lep_e = caja->setup->param->electron;
    
    // neutrino indexes
    int neu_tau = numneu*lep_tau+lep_tau;
    int neu_muon = numneu*lep_muon+lep_muon;
    int neu_e = numneu*lep_e+lep_e;
    
    // branching
    //double br_muon = 0.1736;
    //double br_e = 0.1785;
    // already included in dNdENuTau_Lep
    double br_muon = 1.0;
    double br_e = 1.0;
    
    // convert to complex state
    // neutrino state
    gsl_vector_complex* state = gsl_vector_complex_calloc(caja->numeqn/2);
    state->data = (double*) state_dbl;
    // antineutrino
    gsl_vector_complex* state_aux = gsl_vector_complex_calloc(caja_aux->numeqn/2);
    state_aux->data = (double*) state_dbl_aux;
    
    // load the cross sections and decay spectra
    // for neutrino
    gsl_matrix_complex* dNdENuTau_All = caja->interactions->dNdENuTau_All;
    gsl_matrix_complex* dNdENuTau_Lep = caja->interactions->dNdENuTau_Lep;
    // for antineutrino
    gsl_matrix_complex* dNdENuTau_All_aux = caja_aux->interactions->dNdENuTau_All;
    gsl_matrix_complex* dNdENuTau_Lep_aux = caja_aux->interactions->dNdENuTau_Lep;
        
    for(int e1 = 0; e1 < e_size; e1++){
        gsl_complex reinj_tau_neu_flux_all = gsl_complex_rect(0.0,0.0);
        gsl_complex reinj_tau_neu_flux_lep = gsl_complex_rect(0.0,0.0);
        gsl_complex reinj_tau_aneu_flux_all = gsl_complex_rect(0.0,0.0);
        gsl_complex reinj_tau_aneu_flux_lep = gsl_complex_rect(0.0,0.0);
        // calculate the integral 
        for(int e2 = e1 + 1 ; e2 < e_size; e2++){
            #ifdef TauRegManual_DEBUG
            cout << e1 << " " << e2 << " "<< gsl_matrix_complex_get(dNdENuTau_All,e2,e1).dat[0] << endl;
            #endif
            
            // all decays
                    gsl_complex nutau_integral_term_all = gsl_complex_mul_real(
                                                    gsl_complex_mul(
                                                        gsl_matrix_complex_get(dNdENuTau_All,e2,e1),
                                                        gsl_vector_complex_get(state,SQR(numneu)*e_size + numneu*e2 + lep_tau)
                                                        ),
                                                    (E_range[e2]-E_range[e2-1])
                                                    );
                    
                    gsl_complex nutau_integral_term_lep = gsl_complex_mul_real(
                                                    gsl_complex_mul(
                                                        gsl_matrix_complex_get(dNdENuTau_Lep,e2,e1),
                                                        gsl_vector_complex_get(state,SQR(numneu)*e_size + numneu*e2 + lep_tau)
                                                        ),
                                                    (E_range[e2]-E_range[e2-1])
                                                    );
                    
                    gsl_complex anutau_integral_term_all = gsl_complex_mul_real(
                                                    gsl_complex_mul(
                                                        gsl_matrix_complex_get(dNdENuTau_All_aux,e2,e1),
                                                        gsl_vector_complex_get(state_aux,SQR(numneu)*e_size + numneu*e2 + lep_tau)
                                                        ),
                                                    (E_range[e2]-E_range[e2-1])
                                                    );
                    
                    gsl_complex anutau_integral_term_lep = gsl_complex_mul_real(
                                                    gsl_complex_mul(
                                                        gsl_matrix_complex_get(dNdENuTau_Lep_aux,e2,e1),
                                                        gsl_vector_complex_get(state_aux,SQR(numneu)*e_size + numneu*e2 + lep_tau)
                                                        ),
                                                    (E_range[e2]-E_range[e2-1])
                                                    );
                    
                    reinj_tau_neu_flux_all = gsl_complex_add(nutau_integral_term_all,reinj_tau_neu_flux_all);
                    reinj_tau_neu_flux_lep = gsl_complex_add(nutau_integral_term_lep,reinj_tau_neu_flux_lep);
                    reinj_tau_aneu_flux_all = gsl_complex_add(anutau_integral_term_all,reinj_tau_aneu_flux_all);
                    reinj_tau_aneu_flux_lep = gsl_complex_add(anutau_integral_term_lep,reinj_tau_aneu_flux_lep);
        }        
        // add to current tau neutrino flux
        #ifdef TauRegManual_DEBUG
            cout << "ReinjFlux" << endl;
            cout << e1 << " " << reinj_tau_neu_flux_all.dat[0] << endl;
        #endif
        
        // neutrinos
        // tau neutrinos from tau leptons
        gsl_complex tau_neu_flux = gsl_vector_complex_get(state,SQR(numneu)*e1+neu_tau);
        tau_neu_flux = gsl_complex_add(tau_neu_flux,reinj_tau_neu_flux_all);
        gsl_vector_complex_set(state,SQR(numneu)*e1+neu_tau,tau_neu_flux);
        
        // muon neutrinos from antitau leptons
        gsl_complex muon_neu_flux = gsl_vector_complex_get(state,SQR(numneu)*e1+neu_muon);
        muon_neu_flux = gsl_complex_add(muon_neu_flux,gsl_complex_mul_real(reinj_tau_aneu_flux_lep,br_muon));
        gsl_vector_complex_set(state,SQR(numneu)*e1+neu_muon,muon_neu_flux);
        
        // electron neutrinos from antitau leptons
        gsl_complex e_neu_flux = gsl_vector_complex_get(state,SQR(numneu)*e1+neu_e);
        e_neu_flux = gsl_complex_add(e_neu_flux,gsl_complex_mul_real(reinj_tau_aneu_flux_lep,br_e));
        gsl_vector_complex_set(state,SQR(numneu)*e1+neu_e,e_neu_flux);
        
        // antineutrinos
        // tau neutrinos from tau leptons
        gsl_complex tau_aneu_flux = gsl_vector_complex_get(state_aux,SQR(numneu)*e1+neu_tau);
        tau_aneu_flux = gsl_complex_add(tau_aneu_flux,reinj_tau_aneu_flux_all);
        gsl_vector_complex_set(state_aux,SQR(numneu)*e1+neu_tau,tau_aneu_flux);
        
        // muon neutrinos from antitau leptons
        gsl_complex muon_aneu_flux = gsl_vector_complex_get(state_aux,SQR(numneu)*e1+neu_muon);
        muon_aneu_flux = gsl_complex_add(muon_aneu_flux,gsl_complex_mul_real(reinj_tau_neu_flux_lep,br_muon));
        gsl_vector_complex_set(state_aux,SQR(numneu)*e1+neu_muon,muon_aneu_flux);
        
        // electron neutrinos from antitau leptons
        gsl_complex e_aneu_flux = gsl_vector_complex_get(state_aux,SQR(numneu)*e1+neu_e);
        e_aneu_flux = gsl_complex_add(e_aneu_flux,gsl_complex_mul_real(reinj_tau_neu_flux_lep,br_e));
        gsl_vector_complex_set(state_aux,SQR(numneu)*e1+neu_e,e_aneu_flux);
    }
    
    // remove all the tau lepton and antilepton flux
    for(int e1 = 0; e1 < e_size; e1++){
        gsl_vector_complex_set(state,SQR(numneu)*e_size + numneu*e1 + lep_tau,gsl_complex_rect(0.0,0.0));
        gsl_vector_complex_set(state_aux,SQR(numneu)*e_size + numneu*e1 + lep_tau,gsl_complex_rect(0.0,0.0));
    }
    
    // write back to double representation state
    for (int i = 0 ; i< caja->numeqn; i++){
        state_dbl[i] = state->data[i];
        state_dbl_aux[i] = state_aux->data[i];  
    }
    
    gsl_vector_complex_free(state);
    gsl_vector_complex_free(state_aux);
    
    return GSL_SUCCESS;
}

gsl_matrix_complex* CalculateRHSMatrix(double x,InterOperator* interactions,OscSetup* setup){
    // determine matrix size
    int numneu = setup->param->numneu;
    int numstate = SQR(numneu)+numneu;
    int e_size = setup->E_range.size();
    
    int matrix_size = numstate*e_size;
    // create matrix
    //gsl_matrix_complex* RHS_matrix = gsl_matrix_complex_calloc(matrix_size,matrix_size);
    // zero matrix
    gsl_matrix_complex * RHS_matrix = setup->RHS_matrix;
    gsl_matrix_complex_set_zero(RHS_matrix);
    // fill matrix with hamiltonian operator flavor diagonal terms
    // calculate the hamiltonian conmutator
    setup->track->x = x;
    setup->param->neutype = setup->neutype;
    
    // calculate and store change basis operators for interaction basis
    vector<gsl_matrix_complex*> Uint_vector;
    if (setup->basis == 2){
        // when on interaction basis create evolution operator for H0
        // for each energy
        for(int e = 0; e < e_size; e++){
            // create interaction basis change matrix
            gsl_matrix_complex* Uint = massExpH0(setup->param,setup->E_range[e],
                                            setup->track,setup->massM2);
            Uint_vector.push_back(Uint);
            // update flavor projectors in interaction basis
            for(int flavor = 0; flavor < numneu; flavor++){
                    gsl_matrix_complex* proj = setup->VectorFlavorProjector[e][flavor];
                    gsl_matrix_complex_memcpy(proj,setup->FlavorProjector[flavor]);
                    gsl_matrix_complex_change_basis_UMUC(Uint,proj);
            }
        }
    }
    
    if (setup->oscillation){
        for(int e = 0; e < e_size; e++){
            gsl_matrix_complex *H_current;
            
            if(setup->basis == 0){
                // use flavor basis
                H_current = flavorH(setup->param,setup->E_range[e],
                                    setup->body,setup->track,setup->flavorM2);
            } else if (setup->basis == 1){
                // use mass basis
                H_current = massH(setup->param,setup->E_range[e],
                                  setup->body,setup->track,setup->massM2,setup->U);
            } else if (setup->basis == 2){
                // use interaction basis (with respect to mass basis)
                            
                H_current = massAcc(setup->param,setup->E_range[e],
                                    setup->body,setup->track,setup->U);
                gsl_matrix_complex_change_basis_UMUC(Uint_vector[e],H_current);
            }
            
            // multiply hamiltonian by -i, look at reference.
            gsl_matrix_complex_scale(H_current, gsl_complex_rect(0.0,-1.0));
            
            #ifdef RHS_Hamiltonian_DEBUG
                cout << "== -i X Hamiltonian at x : " << x << " E : " << setup->E_range[e] << " ==" <<  endl;
                for(int i = 0; i < numneu; i++){
                    for(int j = 0; j < numneu; j++){
                        cout << gsl_matrix_complex_get(H_current,i,j).dat[0] <<
                        "+i" << gsl_matrix_complex_get(H_current,i,j).dat[1] << " ";
                }
                cout << endl;
            }
                
                cout << "==" << endl;
            #endif
            
            #ifdef RHS_Matrix_Detail_DEBUG
                cout << "== BEGIN Hamiltonian-RHS-SubMatrix FOR e : " << e << " ==" << endl;
            #endif
            
            // now we calculate i[H,F] which is base invariant
            for(int i = 0; i < SQR(numneu); i++){
                for(int j = 0; j < SQR(numneu); j++){
                    // convert to other indexes
                    int a = i/numneu;
                    int b = i%numneu;
                    int c = j/numneu;
                    int d = j%numneu;
                    
                    #ifdef RHS_Matrix_Osc_DEBUG
                    cout << i << " " << j << endl;
                    cout << a << " " << b << " " << c << " " << d << endl;
                    cout << "RHS_i " << e*SQR(numneu)+i << " RHS_j " << e*SQR(numneu)+j << endl;
                    #endif
                    // implementing the conmutator
                    gsl_matrix_complex_set(RHS_matrix,e*SQR(numneu)+i,e*SQR(numneu)+j,
                                           gsl_complex_sub(
                                           gsl_complex_mul_real(gsl_matrix_complex_get(H_current,a,c),
                                                                (double) KRONECKER(b,d)),
                                           gsl_complex_mul_real(gsl_matrix_complex_get(H_current,d,b),
                                                                (double) KRONECKER(a,c))
                                                          )
                                          );
                }
            }
            
            #ifdef RHS_Matrix_Osc_DEBUG
                cout << "==" << endl;
            #endif
            gsl_matrix_complex_free(H_current);
        }
        
        #ifdef RHS_Matrix_DEBUG
            cout << "== RHS at x : " << x << " ==" <<  endl;
            for(int i = 0; i < matrix_size; i++){
                for(int j = 0; j < matrix_size; j++){
                    cout << gsl_matrix_complex_get(RHS_matrix,i,j).dat[0] <<
                    "+i" << gsl_matrix_complex_get(RHS_matrix,i,j).dat[1] << " ";
                }
                cout << endl;
            }
            cout << "==" << endl;
        #endif
    }

    if (setup->attenuation or setup->nc_inter or setup->cc_inter or setup->tau_regeneration){
        #ifdef RHS_Matrix_IntVector_DEBUG
        if (x > 10.0*setup->param->km ){
            cout << "== RHS Interaction Vector Before Update ==" << endl;
            gsl_vector_complex* Lint_before = interactions->Lint[0];
            for(int e = 0; e < e_size; e++){
                cout << gsl_vector_complex_get(Lint_before,e).dat[0] << " ";
            }
            cout << "=="<< endl;
        }
        #endif
        UpdateInter(interactions,setup);
        #ifdef RHS_Matrix_IntVector_DEBUG
        cout << "Position " << x/setup->param->km << endl;
        cout << "Density " << setup->body->density(setup->track) << endl;
        
        if (x > 10.0*setup->param->km ){
            cout << "== RHS Interaction Vector After ==" << endl;
            gsl_vector_complex* Lint_after = interactions->Lint[0];
            for(int e = 0; e < e_size; e++){
                cout << gsl_vector_complex_get(Lint_after,e).dat[0] << " ";
            }
            cout << "=="<< endl;
        }
        #endif
    }

    // fill matrix with interaction cross terms
    if (setup->attenuation){
        for(int flavor = 0; flavor < MIN(numneu,3); flavor++){
            #ifdef AttenuationOnlyCC
                gsl_vector_complex* Lint = interactions->LCC[flavor];
            #else
                gsl_vector_complex* Lint = interactions->Lint[flavor];
            #endif

            gsl_matrix_complex* RHS_matrix_TMP = gsl_matrix_complex_calloc(matrix_size,matrix_size);
            
            for(int e = 0; e < e_size; e++){
                // getting attenuation matrix
                // get projector into flavor state
                gsl_matrix_complex* attenuation = gsl_matrix_complex_alloc(numneu,numneu);
                
                if(setup->basis == 0 or setup->basis == 1){
                    // get projector
                    gsl_matrix_complex_memcpy(attenuation,setup->FlavorProjector[flavor]);
                } else if(setup->basis == 2){
                    // get projector in interaction basis
                    gsl_matrix_complex_memcpy(attenuation,setup->VectorFlavorProjector[e][flavor]);
                }
                
                // multiply by constant
                gsl_matrix_complex_scale(attenuation,
                                        gsl_complex_mul_real(
                                            gsl_complex_inverse(
                                                gsl_vector_complex_get(Lint,e)),
                                            -0.5)
                                        );

                for(int i = 0; i < SQR(numneu); i++){
                    for(int j = 0; j < SQR(numneu); j++){
                        int a = i/numneu;
                        int b = i%numneu;
                        int c = j/numneu;
                        int d = j%numneu;
                        
                        // implementing the anticonmutator
                        gsl_matrix_complex_set(RHS_matrix_TMP,e*SQR(numneu)+i,e*SQR(numneu)+j,
                                               gsl_complex_add(
                                               gsl_complex_mul_real(gsl_matrix_complex_get(attenuation,a,c),
                                                                    (double) KRONECKER(b,d)),
                                               gsl_complex_mul_real(gsl_matrix_complex_get(attenuation,d,b),
                                                                    (double) KRONECKER(a,c))
                                                              )
                                              );
                    }
                }
                #ifdef RHS_Matrix_Int_DEBUG
                cout << "== RHS Int Flavor " << flavor << " ==" << endl;
                for(int i = 0; i < matrix_size; i++){
                    for(int j = 0; j < matrix_size; j++){
                        cout << gsl_matrix_complex_get(RHS_matrix_TMP,i,j).dat[0] <<
                        "+i" << gsl_matrix_complex_get(RHS_matrix_TMP,i,j).dat[1] << " ";
                    }
                    cout << endl;
                }
                cout << "==" << endl;
                #endif
                gsl_matrix_complex_free(attenuation);
            }
            
            #ifdef CUDA_ENABLE
                CUDAComplexMatrixAdd(RHS_matrix,RHS_matrix_TMP);
            #else
                gsl_matrix_complex_add(RHS_matrix,RHS_matrix_TMP);
            #endif
            gsl_matrix_complex_free(RHS_matrix_TMP);
                        
        }
        #ifdef RHS_Matrix_DEBUG
            cout << "== RHS at x : " << x << " ==" <<  endl;
            for(int i = 0; i < matrix_size; i++){
                for(int j = 0; j < matrix_size; j++){
                    cout << gsl_matrix_complex_get(RHS_matrix,i,j).dat[0] <<
                    "+i" << gsl_matrix_complex_get(RHS_matrix,i,j).dat[1] << " ";
                }
                cout << endl;
            }
            cout << "==" << endl;
        #endif
    }

    if (setup->nc_inter) {
        gsl_vector_complex* sigma_NC = interactions->sigma_NC[0];
        gsl_vector_complex* LNC = interactions->LNC[0];
        gsl_matrix_complex* dsignudE_NC = interactions->dsignudE_NC[0];
        
        gsl_matrix_complex* RHS_matrix_TMP = gsl_matrix_complex_calloc(matrix_size,matrix_size);
        
        for(int e1 = 0; e1 < e_size; e1++){
            for(int e2 = e1 + 1 ; e2 < e_size; e2++){
                gsl_complex NC_integral_term = gsl_complex_mul_real(
                                                    gsl_complex_div(gsl_matrix_complex_get(dsignudE_NC,e2,e1),
                                                        gsl_complex_mul(
                                                            gsl_vector_complex_get(LNC,e2),
                                                            gsl_vector_complex_get(sigma_NC,e2)
                                                                       )
                                                                   ),
                                                (setup->E_range[e2]-setup->E_range[e2-1])
                                                );
                
                #ifdef RHS_Matrix_NC_DEBUG
                    cout << "e1 " << e1 << " e2 " << e2 << endl;
                    cout << NC_integral_term.dat[0] << "+i" << NC_integral_term.dat[1] << endl;
                #endif
                
                // this term is flavor blind
                for(int i = 0; i < SQR(numneu); i++){
                    gsl_matrix_complex_set(RHS_matrix_TMP,
                                           e1*SQR(numneu)+i,e2*SQR(numneu)+i,
                                           NC_integral_term);       
                }
            }
        }
        
        #ifdef CUDA_ENABLE
            CUDAComplexMatrixAdd(RHS_matrix,RHS_matrix_TMP);
        #else
            gsl_matrix_complex_add(RHS_matrix,RHS_matrix_TMP);
        #endif
        
        gsl_matrix_complex_free(RHS_matrix_TMP);
        
        #ifdef RHS_Matrix_NC_DEBUG
            cout << "== RHS at x : " << x << " ==" <<  endl;
            for(int i = 0; i < matrix_size; i++){
                for(int j = 0; j < matrix_size; j++){
                    cout << gsl_matrix_complex_get(RHS_matrix,i,j).dat[0] <<
                    "+i" << gsl_matrix_complex_get(RHS_matrix,i,j).dat[1] << " ";
                }
                cout << endl;
            }
            cout << "==" << endl;
        #endif
    }
    
    if (setup->cc_inter) {
        // nothing
    }
    
    if(setup->muon){
        int lep_muon = setup->param->muon;
        int neu_muon = numneu*lep_muon+lep_muon;
        
        // define some pointers
        gsl_vector_complex* sigma_muon_CC = interactions->sigma_CC[lep_muon];
        gsl_matrix_complex* dsignudE_muon_CC = interactions->dsignudE_CC[lep_muon];
        
        gsl_vector_complex* LCC_muon = interactions->LCC[lep_muon];
        
        ///////////////////////////// lepton term ///////////////////////////
        // define temporary mixing matrix
        gsl_matrix_complex* RHS_matrix_lep_TMP = gsl_matrix_complex_calloc(matrix_size,matrix_size);

        // integral term
        for(int e1 = 0; e1 < e_size; e1++){
            for(int e2 = e1 + 1 ; e2 < e_size; e2++){
                gsl_complex muon_CC_integral_term = gsl_complex_mul_real(
                                                    gsl_complex_div(gsl_matrix_complex_get(dsignudE_muon_CC,e2,e1),
                                                        gsl_complex_mul(
                                                            gsl_vector_complex_get(LCC_muon,e2),
                                                            gsl_vector_complex_get(sigma_muon_CC,e2)
                                                                       )
                                                                   ),
                                                (setup->E_range[e2]-setup->E_range[e2-1])
                                                );
                
                if(setup->basis == 0){
                    // selecting the muon component makes the TraceProyect
                    gsl_matrix_complex_set(RHS_matrix_lep_TMP,
                                            e_size*SQR(numneu)+numneu*e1+lep_muon,
                                            e2*SQR(numneu)+neu_muon, 
                                            muon_CC_integral_term);
                } else if (setup->basis == 1){
                    gsl_matrix_complex* muon_term_matrix = gsl_matrix_complex_alloc(numneu,numneu);
                    gsl_matrix_complex_memcpy(muon_term_matrix,setup->FlavorProjector[lep_muon]);
                    gsl_matrix_complex_scale(muon_term_matrix,muon_CC_integral_term);
                    
                    for(int i = 0; i < SQR(numneu); i++){
                        gsl_matrix_complex_set(RHS_matrix_lep_TMP,
                                                e_size*SQR(numneu)+numneu*e1+lep_muon,
                                                e2*SQR(numneu)+i, 
                                                gsl_matrix_complex_get(muon_term_matrix,i/numneu,i%numneu)
                                              );
                    }
                    
                    gsl_matrix_complex_free(muon_term_matrix);
                } else if (setup->basis == 2){
                    gsl_matrix_complex* muon_term_matrix = gsl_matrix_complex_alloc(numneu,numneu);
                    gsl_matrix_complex_memcpy(muon_term_matrix,setup->VectorFlavorProjector[e1][lep_muon]);
                    gsl_matrix_complex_scale(muon_term_matrix,muon_CC_integral_term);
                    
                    for(int i = 0; i < SQR(numneu); i++){
                        gsl_matrix_complex_set(RHS_matrix_lep_TMP,
                                                e_size*SQR(numneu)+numneu*e1+lep_muon,
                                                e2*SQR(numneu)+i, 
                                                gsl_matrix_complex_get(muon_term_matrix,i/numneu,i%numneu)
                                              );
                    }
                    
                    gsl_matrix_complex_free(muon_term_matrix);
                }
            }
        }
        
        // adding to main RHS matrix
        #ifdef CUDA_ENABLE
            CUDAComplexMatrixAdd(RHS_matrix,RHS_matrix_lep_TMP);
        #else
            gsl_matrix_complex_add(RHS_matrix,RHS_matrix_lep_TMP);
        #endif
        gsl_matrix_complex_free(RHS_matrix_lep_TMP);
    }
    
    if (setup->tau_regeneration and numneu > setup->param->tau) {
        int lep_tau = setup->param->tau;
        int neu_tau = numneu*lep_tau+lep_tau;
        
        // define some pointers
        gsl_vector_complex* sigma_tau_CC = interactions->sigma_CC[lep_tau];
        gsl_matrix_complex* dsignudE_tau_CC = interactions->dsignudE_CC[lep_tau];
        
        gsl_vector_complex* LCC_tau = interactions->LCC[lep_tau];
        gsl_vector_complex* LTauDecay = interactions->LTauDecay;
        
        gsl_matrix_complex* dNdENuTau_All = interactions->dNdENuTau_All;
        //gsl_matrix_complex* dNdENuTau_Lep = interactions->dNdENuTau_Lep;
        
        ///////////////////////////// lepton term ///////////////////////////
        // define temporary mixing matrix
        gsl_matrix_complex* RHS_matrix_lep_TMP = gsl_matrix_complex_calloc(matrix_size,matrix_size);

        // integral term
        for(int e1 = 0; e1 < e_size; e1++){
            for(int e2 = e1 + 1 ; e2 < e_size; e2++){
                
                #ifdef RHS_Matrix_Tau_DEBUG
                    cout << "Tau term" << endl;
                    cout << e2 << " " << gsl_vector_complex_get(LCC_tau,e2).dat[0] << endl;
                    cout << e1 << " " <<  e2 << " " << gsl_complex_div(gsl_matrix_complex_get(dsignudE_tau_CC,e2,e1),
                                                                    gsl_vector_complex_get(sigma_tau_CC,e2)).dat[0] << endl;
                #endif
                
                gsl_complex tau_CC_integral_term = gsl_complex_mul_real(
                                                    gsl_complex_div(gsl_matrix_complex_get(dsignudE_tau_CC,e2,e1),
                                                        gsl_complex_mul(
                                                            gsl_vector_complex_get(LCC_tau,e2),
                                                            gsl_vector_complex_get(sigma_tau_CC,e2)
                                                                       )
                                                                   ),
                                                (setup->E_range[e2]-setup->E_range[e2-1])
                                                );
                // check this
                if(setup->basis == 0){
                    // selecting the tau component makes the TraceProyect
                    gsl_matrix_complex_set(RHS_matrix_lep_TMP,
                                            e_size*SQR(numneu)+numneu*e1+lep_tau,
                                            e2*SQR(numneu)+neu_tau, 
                                            tau_CC_integral_term);
                } else if (setup->basis == 1){
                    gsl_matrix_complex* tau_term_matrix = gsl_matrix_complex_alloc(numneu,numneu);
                    gsl_matrix_complex_memcpy(tau_term_matrix,setup->FlavorProjector[lep_tau]);
                    gsl_matrix_complex_scale(tau_term_matrix,tau_CC_integral_term);
                    
                    for(int i = 0; i < SQR(numneu); i++){
                        gsl_matrix_complex_set(RHS_matrix_lep_TMP,
                                                e_size*SQR(numneu)+numneu*e1+lep_tau,
                                                e2*SQR(numneu)+i, 
                                                gsl_matrix_complex_get(tau_term_matrix,i/numneu,i%numneu)
                                              );
                    }
                    
                    gsl_matrix_complex_free(tau_term_matrix);
                } else if (setup->basis == 2){
                    gsl_matrix_complex* tau_term_matrix = gsl_matrix_complex_alloc(numneu,numneu);
                    gsl_matrix_complex_memcpy(tau_term_matrix,setup->VectorFlavorProjector[e1][lep_tau]);
                    gsl_matrix_complex_scale(tau_term_matrix,tau_CC_integral_term);
                    
                    for(int i = 0; i < SQR(numneu); i++){
                        gsl_matrix_complex_set(RHS_matrix_lep_TMP,
                                                e_size*SQR(numneu)+numneu*e1+lep_tau,
                                                e2*SQR(numneu)+i, 
                                                gsl_matrix_complex_get(tau_term_matrix,i/numneu,i%numneu)
                                              );
                    }
                    
                    gsl_matrix_complex_free(tau_term_matrix);
                }
                
            }
        }
        
        #ifdef RHS_Matrix_Unitarity_DEBUG
            for(int e1 = 0; e1 < e_size; e1++){
                gsl_complex tau_CC_integral = gsl_complex_rect(0.0,0.0);
                for(int e2 = 0; e2 < e1; e2++){
                    gsl_complex tau_CC_integral_term = gsl_complex_mul_real(
                                                        gsl_complex_div(gsl_matrix_complex_get(dsignudE_tau_CC,e1,e2),
                                                            gsl_complex_mul(
                                                                gsl_vector_complex_get(LCC_tau,e1),
                                                                gsl_vector_complex_get(sigma_tau_CC,e1)
                                                                           )
                                                                       ),
                                                    (setup->E_range[e2]-setup->E_range[e2-1])
                                                    );
                    
                    tau_CC_integral = gsl_complex_add(tau_CC_integral,tau_CC_integral_term);
                }
                
                gsl_complex max_error = gsl_complex_mul_real(
                                                        gsl_complex_div(gsl_matrix_complex_get(dsignudE_tau_CC,e1,0),
                                                            gsl_complex_mul(
                                                                gsl_vector_complex_get(LCC_tau,e1),
                                                                gsl_vector_complex_get(sigma_tau_CC,e1)
                                                                           )
                                                                       ),
                                                    (setup->E_range[0])
                                                    );
                
                double error = max_error.dat[0];
                double width_1 = 1.0/gsl_vector_complex_get(LCC_tau,e1).dat[0];
                double width_2 = tau_CC_integral.dat[0];
                
                cout << "E : " << setup->E_range[e1]*1.0e-9 << " " << width_1 << " " << width_2 << " " << width_2/width_1  << " " << width_2-width_1 << " " << error << endl;
            }
        #endif
        
        #ifndef ManualTauReinjection
            // decay rate term        
            for(int e = 0; e < e_size; e++){
                int entry = e_size*SQR(numneu)+numneu*e+lep_tau;
                gsl_matrix_complex_set(RHS_matrix_lep_TMP,entry,entry,
                                        gsl_complex_div(
                                                        gsl_complex_rect(-1.0,0.0),
                                                        gsl_vector_complex_get(LTauDecay,e))
                                       );
            }
        #endif
        
        #ifdef RHS_Matrix_tau_lepton_DEBUG
            cout << "== RHS at x : " << x << " ==" <<  endl;
            for(int i = 0; i < matrix_size; i++){
                for(int j = 0; j < matrix_size; j++){
                    cout << gsl_matrix_complex_get(RHS_matrix_lep_TMP,i,j).dat[0] <<
                    "+i" << gsl_matrix_complex_get(RHS_matrix_lep_TMP,i,j).dat[1] << " ";
                }
                cout << endl;
            }
            cout << "==" << endl;
        #endif
        
        // adding to main RHS matrix
        #ifdef CUDA_ENABLE
            CUDAComplexMatrixAdd(RHS_matrix,RHS_matrix_lep_TMP);
        #else
            gsl_matrix_complex_add(RHS_matrix,RHS_matrix_lep_TMP);
        #endif
        gsl_matrix_complex_free(RHS_matrix_lep_TMP);
        
        ///////////////////////////// neutrino term ///////////////////////////
        gsl_matrix_complex* RHS_matrix_neu_TMP = gsl_matrix_complex_calloc(matrix_size,matrix_size);        
        
        // integral term
        #ifndef ManualTauReinjection
            for(int e1 = 0; e1 < e_size; e1++){
                for(int e2 = e1 + 1 ; e2 < e_size; e2++){
                    gsl_complex nutau_integral_term = gsl_complex_mul_real(
                                                        gsl_complex_div(gsl_matrix_complex_get(dNdENuTau_All,e2,e1),
                                                                        gsl_vector_complex_get(LTauDecay,e2)
                                                                        ),
                                                    (setup->E_range[e2]-setup->E_range[e2-1])
                                                    );
                    if(setup->basis == 0){
                        gsl_matrix_complex_set(RHS_matrix_neu_TMP,
                                               e1*SQR(numneu)+neu_tau,
                                               e_size*SQR(numneu)+numneu*e2+lep_tau,
                                               nutau_integral_term);
                    } else if (setup->basis == 1){
                    // need to work here     
                    } else if (setup->basis == 2){
                    // need to work here
                    }
                }
            }
        #endif
        
        #ifdef RHS_Matrix_tau_neutrino_DEBUG
            cout << "== RHS at x : " << x << " ==" <<  endl;
            for(int i = 0; i < matrix_size; i++){
                for(int j = 0; j < matrix_size; j++){
                    cout << gsl_matrix_complex_get(RHS_matrix_neu_TMP,i,j).dat[0] <<
                    "+i" << gsl_matrix_complex_get(RHS_matrix_neu_TMP,i,j).dat[1] << " ";
                }
                cout << endl;
            }
            cout << "==" << endl;
        #endif
        
        
        // adding to main RHS matrix
        #ifdef CUDA_ENABLE
            CUDAComplexMatrixAdd(RHS_matrix,RHS_matrix_neu_TMP);
        #else
            gsl_matrix_complex_add(RHS_matrix,RHS_matrix_neu_TMP);
        #endif
        gsl_matrix_complex_free(RHS_matrix_neu_TMP);
    }
    
    if (setup->basis == 2){
        for(int e = 0; e < e_size; e++){
            gsl_matrix_complex_free(Uint_vector[e]);
        }
    }
    
    return RHS_matrix;
}

int RHS_RHO(double x ,const double *state_dbl_in,double *state_dbl_out,void *box){
    // Format container
    Box* caja = (Box*) box;
    // Initialize state vectors
    gsl_vector_complex* state_in = gsl_vector_complex_calloc(caja->numeqn/2);
    state_in->data = (double*) state_dbl_in;

    #ifdef RHS_RHO_DEBUG
        cout << "==RHS IN state==" << endl;
        /*
        for (int i = 0 ; i< caja->numeqn; i++){
            cout << state_dbl_in[i] << " ";     
        }
        cout << endl;
        */
        cout << "state size "<< state_in->size << endl;
        cout << "state value "<< endl;
        gsl_vector_complex_fprintf(stdout,state_in,"%g");
        cout << "==" << endl;
    #endif
    
    // Create complex vector out 
    gsl_vector_complex* state_out = gsl_vector_complex_calloc(state_in->size);
    // Calculate RHS matrix
    CalculateRHSMatrix(x,caja->interactions,caja->setup);
    //gsl_matrix_complex* RHS_matrix = CalculateRHSMatrix(x,caja->interactions,caja->setup);
    #ifdef RHS_RHO_DEBUG
        cout << "== RHS Matrix : x = " << x/caja->setup->param->km << " ==" <<  endl;
        for(unsigned int i = 0; i < caja->setup->RHS_matrix->size1; i++){
            for(unsigned int j = 0; j < caja-.setup->RHS_matrix->size1; j++){
                cout << gsl_matrix_complex_get(caja->setup->RHS_matrix,i,j).dat[0] <<
                "+i" << gsl_matrix_complex_get(caja->setup->RHS_matrix,i,j).dat[1] << " ";
            }
            cout << endl;
        }
        cout << "==" << endl;
    #endif
    // Multiply by state vector -> GPU
    gsl_complex one = gsl_complex_rect(1.0,0.0);
    gsl_complex zero = gsl_complex_rect(0.0,0.0);
    gsl_blas_zgemv(CblasNoTrans,one,caja->setup->RHS_matrix,state_in,zero,state_out);
    #ifdef RHS_RHO_DEBUG
        cout << "==RHS OUT state==" << endl;
        cout << "state size "<< state_out->size << endl;
        cout << "state value "<< endl;
        gsl_vector_complex_fprintf(stdout,state_out,"%g");
        cout << "==" << endl;
    #endif
    // lets go back to the double representation
    //state_dbl_out = state_out->data;
    for (int i = 0 ; i< caja->numeqn; i++){
            state_dbl_out[i] = state_out->data[i];     
    }
    
    #ifdef RHS_RHO_DEBUG
    /*
        for (int i = 0 ; i< caja->numeqn; i++){
            cout << state_dbl_out[i] << " ";     
        }
        cout << endl;
    */
    #endif
    
    // free memory
    gsl_vector_complex_free(state_in);
    gsl_vector_complex_free(state_out);
    //gsl_matrix_complex_free(RHS_matrix);
    
    return GSL_SUCCESS;
}

int JAC_RHO(double x ,const double *state_dbl_in,double *dfdx, double dfdt[],void *box){
    // Format container
    Box* caja = (Box*) box;
    // calculate jacobian as gsl matrix complex
    gsl_matrix_complex* JAC_matrix = CalculateRHSMatrix(x,caja->interactions,caja->setup);
    // pass jacobian to double array
    int numeqn = caja->numeqn;
    for(int i = 0; i < numeqn; i++){
        for(int j = 0; j < numeqn; j = j + 2){
            dfdx[i*numeqn + j + 0]   = gsl_matrix_complex_get(JAC_matrix,i,j).dat[0];
            dfdx[i*numeqn + j + 1] = gsl_matrix_complex_get(JAC_matrix,i,j).dat[1];
        }
    }    
    // free gsl complex matrix
    gsl_matrix_complex_free(JAC_matrix);
    
    // set all time derivative to zero
    for(int i = 0; i < numeqn; i++){
        dfdt[i] = 0.0;
    }
    
    return GSL_SUCCESS;
}

void ConfigureSetup(OscSetup* setup){
    int numneu = setup->param->numneu;
    int numstate = SQR(numneu)+numneu;
    int e_size = setup->E_range.size();
    
    int matrix_size = numstate*e_size;
    
    // create setup matrices
    //gsl_matrix_complex_memcpy(setup->U,MixMatrix(setup->param));
    setup->U = MixMatrix(setup->param);
    setup->flavorM2 = flavorM2(setup->param);
    setup->massM2 = massM2(setup->param);
    // projector
    setup->FlavorProjector = InitProjectorVector(setup);
    // intermediate operation matrices
    setup->RHS_matrix = gsl_matrix_complex_alloc(matrix_size,matrix_size);
    
    if (setup->basis == 2){
        vector< vector<gsl_matrix_complex*> > VectorFlavorProjector;
        for(unsigned int e = 0; e < setup->E_range.size(); e++){
            vector<gsl_matrix_complex*> FlavorProjectors;
            for(int flavor = 0; flavor < numneu; flavor++){
                FlavorProjectors.push_back(gsl_matrix_complex_alloc(numneu,numneu));
            }
            VectorFlavorProjector.push_back(FlavorProjectors);
        }
        setup->VectorFlavorProjector = VectorFlavorProjector;
    }
}

VectorFluxState** CalNeuOscInt(VectorFluxState* state_vector_in[2], OscSetup* setup,double abs_error,double rel_error,bool optimization = true){
    int gsl_status = GSL_SUCCESS;
       
    int numneu = setup->param->numneu;
    int e_size = setup->E_range.size();
    int numeqn = 2*e_size*(numneu + SQR(numneu));
    
    bool neu_and_aneu = false;
    
    if( setup->neutype == 0 or setup->neutype == 1) {
        // one and only one neutrino type
        setup->param->neutype = setup->neutype;
        neu_and_aneu = false;
    } else {
        // neutrinos and antineutrinos
        neu_and_aneu = true;
    }
    
    // initial setup. containter setup
    Box* caja = new Box;
    caja->setup = setup;
    if(neu_and_aneu){
        setup->neutype = 0;
        setup->param->neutype = 0;
    }
    ConfigureSetup(setup);
    caja->interactions = InitInter(setup);
    caja->numeqn = numeqn;
    
    // aux initial setup
    Box* caja_aux = new Box;
    // clone setup object
    OscSetup* setup_aux = new OscSetup(*setup);
    caja_aux->setup = setup_aux;
    if(neu_and_aneu){
        setup_aux->neutype = 1;
        setup_aux->param->neutype = 1;
    }
    ConfigureSetup(setup_aux);
    caja_aux->interactions = InitInter(setup_aux);
    caja_aux->numeqn = numeqn;
    
    // setting up GSL ODE solver
    gsl_odeiv2_step *s = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45, numeqn);
    //gsl_odeiv2_step *s = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rk8pd, numeqn);
    gsl_odeiv2_control *c = gsl_odeiv2_control_y_new(rel_error,abs_error);
    gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc(numeqn);
    
    // ODE system
    gsl_odeiv2_system sys = {&RHS_RHO, &JAC_RHO, numeqn, caja};
    gsl_odeiv2_system sys_aux = {&RHS_RHO, &JAC_RHO, numeqn, caja_aux};
    //gsl_odeiv2_system sys = {&RHS_RHO, NULL, numeqn, caja};
    
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
        gsl_odeiv2_step_set_driver(s,d_aux);
        //gsl_odeiv2_driver_set_hmin(d_aux,h_min);
        gsl_odeiv2_driver_set_hmax(d_aux,h_max);
    #endif  
    
    #ifdef CalNeuOscInt_DEBUG
        printf("GSL paramers :\n");
        printf("x_ini : %g [km]\n", x_ini/km);
        printf("x_end : %g [km]\n", x_end/km);
        printf("h : %g [km]\n", h/km);
        printf("h_min : %g [km]\n", h_min/km);      
    #endif
    
    // initial position
    x = x_ini;
    x_aux = x_ini;

    // setup initial state
    // convert statevector to complex vector
    gsl_vector_complex* state_in = ConvertVectorStateToVector(state_vector_in[0],caja->setup);

    gsl_vector_complex* state_in_aux = NULL;
    if (neu_and_aneu){
        state_in_aux = ConvertVectorStateToVector(state_vector_in[1],caja_aux->setup);
    }

    // convert to double pointer
    //double* state_dbl = state_in->data;
    double state_dbl[numeqn];
    double state_dbl_aux[numeqn];
    for(int i=0;i<numeqn;i++){
        state_dbl[i]=state_in->data[i];
        if(neu_and_aneu){
            state_dbl_aux[i]=state_in_aux->data[i];
        } else {
            state_dbl_aux[i]=state_in->data[i];  
       }
    }
    
    #ifdef CalNeuOscInt_DEBUG
        int count = 0;
        int count_step = 10;

    #endif
    
    #ifdef USE_FIX_STEP
        h = 1.0*km;
    #endif

    while (x < x_end){
        #ifdef RK_Unitarity_DEBUG
            // add all
            double total_flux = 0.0;
            for(int e1 = 1; e1 < e_size; e1++){
                int flv = setup->param->tau;
                //for(int flv = 0; flv < numneu; flv++){                    
                    total_flux += state_dbl[e1*SQR(numneu)*2+(flv+flv*numneu)*2]*(setup->E_range[e1]-setup->E_range[e1-1]);
                    total_flux += state_dbl[e_size*SQR(numneu)*2+e1*numneu*2+(flv)*2]*(setup->E_range[e1]-setup->E_range[e1-1]);
                //}
            }
            cout << "TFlux " << total_flux << endl;
        #endif
        
        #ifdef USE_GSL_ODE_DRIVER
            double x_inter = x + 100.0*km;
            gsl_status = gsl_odeiv2_driver_apply(d,&x,x_inter,state_dbl);
            if( neu_and_aneu ){
                gsl_status = gsl_odeiv2_driver_apply(d_aux,&x_aux,x_inter,state_dbl_aux);
            }
        #elseif USE_FIX_STEP
            gsl_status = gsl_odeiv2_evolve_apply_fixed_step(e,c,s,&sys,&x,h,state_dbl);
            if( neu_and_aneu ){
                gsl_status = gsl_odeiv2_evolve_apply_fixed_step(e,c,s,&sys_aux,&x_aux,h,state_dbl_aux);
            }
        #else
            gsl_status = gsl_odeiv2_evolve_apply(e,c,s,&sys,&x,x_end,&h,state_dbl);
            if( neu_and_aneu ){
                cout << "aqui" << endl;
                gsl_status = gsl_odeiv2_evolve_apply(e,c,s,&sys_aux,&x_aux,x_end,&h,state_dbl_aux);
            }
        #endif

        #ifdef ManualTauReinjection
            if(setup->tau_regeneration and not neu_and_aneu){
                TauLeptonReinjection(state_dbl,caja);
            } else {
                TauLeptonReinjection(state_dbl,caja,state_dbl_aux,caja_aux);
            }
        #endif
        
        #ifdef CalNeuOscInt_DEBUG
            if(count%count_step == 0){
                printf("x_current : %g %g [km]\n", x/km,x_aux/km);
                //printf("step : %g [km]\n", h/km);
            }
        #endif
        
        #ifdef CalNeuOscIntStep_DEBUG
            if(count%count_step == 0 ) {
                cout << "state_current : " << endl;
                for (int i = 0 ; i< caja->numeqn; i++){
                    cout << state_dbl[i] << " ";     
                }
                cout << endl;
            }
        #endif

        if( gsl_status != GSL_SUCCESS ){
            fprintf(stderr,"CalNeuOscGSL: Error in GSL ODE solver,\n");
            break;
        }
        
        if(h < h_min){
            h = h_min;
        }
        
        if(h > h_max){
            h = h_max;
        }
        
        #ifdef CalNeuOscInt_DEBUG
            count++;
        #endif
    }
    
    if (e) { gsl_odeiv2_evolve_free(e);      e = NULL; }
    if (c) { gsl_odeiv2_control_free(c);     c = NULL; }
    if (s) { gsl_odeiv2_step_free(s);        s = NULL; }
    #ifdef USE_GSL_ODE_DRIVER
    if (d) { gsl_odeiv2_driver_free(d);        d = NULL; }
    if (d_aux) { gsl_odeiv2_driver_free(d_aux);        d_aux = NULL; }    
    #endif

    // convert back to complex vector
    gsl_vector_complex* state_out = gsl_vector_complex_calloc(caja->numeqn/2);
    state_out->data = (double*) state_dbl;
    
    gsl_vector_complex* state_out_aux = gsl_vector_complex_calloc(caja->numeqn/2);
    state_out_aux->data = (double*) state_dbl_aux;
    // convert to vector state
    VectorFluxState *state_vector_out[2];// = new VectorFluxState[2];
    
    state_vector_out[0] = ConvertVectorToVectorState(state_out,caja->setup);
    if(neu_and_aneu){
        state_vector_out[1] = ConvertVectorToVectorState(state_out_aux,caja_aux->setup);
    }    

    return state_vector_out;    
}
