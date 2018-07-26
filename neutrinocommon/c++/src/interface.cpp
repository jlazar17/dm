#include "interface.h"

int InitLoadFlux(VectorFluxState* state[2],string path,OscSetup* setup){
    VectorFluxState* state_neu = new VectorFluxState;
    VectorFluxState* state_aneu = new VectorFluxState;
    
    // create energy structure
    vector<double> E_range;

    // load data
    Table dat_array = quickread(path);
    int e_size = dat_array.size();
    int numneu = (dat_array[0].size()-1)/2;
    
    for(int i = 0; i < e_size; i++){
        Row dat = dat_array[i];
        
        // save energy
        E_range.push_back(dat[0]);
        
        // neutrinos
        gsl_matrix_complex* F_nu = gsl_matrix_complex_calloc(numneu,numneu);
        gsl_matrix_complex* F_anu = gsl_matrix_complex_calloc(numneu,numneu);
        for(int j = 0; j < numneu; j++){

        
        gsl_matrix_complex_set(F_nu,
                                   j,j,
                                   gsl_complex_rect(dat[j+1],0));

        gsl_matrix_complex_set(F_anu,
                                   j,j,
                                   gsl_complex_rect(dat[j+numneu+1],0));
        
        }
        
        // leptons
        gsl_matrix_complex* F_le = gsl_matrix_complex_calloc(numneu,1);
        gsl_matrix_complex* F_ale = gsl_matrix_complex_calloc(numneu,1);
        
        // put into flux state objects
        
        FluxState* sneu = new FluxState;
        sneu->F_nu = F_nu;
        sneu->F_le = F_le;
        // putting into state vector
        state_neu->push_back(sneu);
        
        FluxState* saneu = new FluxState;
        saneu->F_nu = F_anu;
        saneu->F_le = F_ale;
        // putting into state vector
        state_aneu->push_back(saneu);
    }
    
    // put into structure
    state[0] = state_neu;
    state[1] = state_aneu;
    
    // reconfigure setup
    setup->param->numneu = numneu;
    // remember this
    //setup->neutype = 2;
    setup->E_range = E_range;
    
    // go gundam
    return GSL_SUCCESS;
}

int InitLoadFlux(VectorFluxState* state[2],Table dat_array,OscSetup* setup){
    VectorFluxState* state_neu = new VectorFluxState;
    VectorFluxState* state_aneu = new VectorFluxState;
    
    // create energy structure
    vector<double> E_range;
    
    // load data
    int e_size = dat_array.size();
    int numneu = (dat_array[0].size()-1)/2;
    
    for(int i = 0; i < e_size; i++){
        Row dat = dat_array[i];
        
        // save energy
        E_range.push_back(dat[0]);
        
        // neutrinos
        gsl_matrix_complex* F_nu = gsl_matrix_complex_calloc(numneu,numneu);
        gsl_matrix_complex* F_anu = gsl_matrix_complex_calloc(numneu,numneu);
        for(int j = 0; j < numneu; j++){
            
            
            gsl_matrix_complex_set(F_nu,
                                   j,j,
                                   gsl_complex_rect(dat[j+1],0));
            
            gsl_matrix_complex_set(F_anu,
                                   j,j,
                                   gsl_complex_rect(dat[j+numneu+1],0));
            
        }
        
        //gsl_matrix_complex_fprintf(stdout,F_nu,"%g");
        //gsl_matrix_complex_fprintf(stdout,F_anu,"%g");
        
        // leptons
        gsl_matrix_complex* F_le = gsl_matrix_complex_calloc(numneu,1);
        gsl_matrix_complex* F_ale = gsl_matrix_complex_calloc(numneu,1);
        
        // put into flux state objects
        
        FluxState* sneu = new FluxState;
        sneu->F_nu = F_nu;
        sneu->F_le = F_le;
        // putting into state vector
        state_neu->push_back(sneu);
        
        FluxState* saneu = new FluxState;
        saneu->F_nu = F_anu;
        saneu->F_le = F_ale;
        // putting into state vector
        state_aneu->push_back(saneu);
    }
    
    // put into structure
    state[0] = state_neu;
    state[1] = state_aneu;
    
    // reconfigure setup
    setup->param->numneu = numneu;
    // remember this
    //setup->neutype = 2;
    setup->E_range = E_range;
    
    // go gundam
    return GSL_SUCCESS;
}

VectorFluxState* InitGeneralFlux(FluxBox* flux,OscSetup* setup,int neutype){
    VectorFluxState* state_vector  = new VectorFluxState;
    
    int e_size = setup->E_range.size();
    // create initial state : neutrino density matrix
    int numneu = setup->param->numneu;
    
    // neutrino indeces
    int nue = setup->param->electron;
    int numu = setup->param->muon;
    int nutau = setup->param->tau;
    
    for(int e=0; e<e_size; e++){
        // initialize neutrino flux density matrix
        gsl_matrix_complex* F_nu = gsl_matrix_complex_calloc(numneu,numneu);
        if (neutype == 0){
            gsl_matrix_complex_set(F_nu,
                                   nue,nue,
                                   gsl_complex_rect(flux->nue_values[e],0));
        
            gsl_matrix_complex_set(F_nu,
                                   numu,numu,
                                   gsl_complex_rect(flux->numu_values[e],0));
                
            gsl_matrix_complex_set(F_nu,
                                   nutau,nutau,
                                   gsl_complex_rect(flux->nutau_values[e],0));
        } else {
            gsl_matrix_complex_set(F_nu,
                                   nue,nue,
                                   gsl_complex_rect(flux->anue_values[e],0));
        
            gsl_matrix_complex_set(F_nu,
                                   numu,numu,
                                   gsl_complex_rect(flux->anumu_values[e],0));
                
            gsl_matrix_complex_set(F_nu,
                                   nutau,nutau,
                                   gsl_complex_rect(flux->anutau_values[e],0));
        }
        
        
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

int SetOscParamMask(double paramlist[],PhysConst* param){
    // set parameters given as a list
    // LIST ORDER IS VERY IMPORTANT
    // number of neutrinos
    param->numneu = paramlist[0];
    // mixing angles
    param->th12 = paramlist[1];
    param->th13 = paramlist[2];
    param->th23 = paramlist[3];
    param->th14 = paramlist[4];
    param->th24 = paramlist[5];
    param->th34 = paramlist[6];
    param->th15 = paramlist[7];
    param->th25 = paramlist[8];
    param->th35 = paramlist[9];
    param->th45 = paramlist[10];
    param->th16 = paramlist[11];
    param->th26 = paramlist[12];
    param->th36 = paramlist[13];
    param->th46 = paramlist[14];
    param->th56 = paramlist[15];
    // square mass differences
    param->dm21sq = paramlist[16];
    param->dm31sq = paramlist[17];
    param->dm41sq = paramlist[18];
    param->dm51sq = paramlist[19];
    param->dm61sq = paramlist[20];
    // cp phases
    param->delta1 = paramlist[21];
    param->delta2 = paramlist[22];
    param->delta3 = paramlist[23];
    // neutrino species
    param->electron = paramlist[24];
    param->muon = paramlist[25];
    param->tau = paramlist[26];
    param->sterile1 = paramlist[27];
    param->sterile2 = paramlist[28];
    param->sterile3 = paramlist[29];
    
    return GSL_SUCCESS;
}

/*
int SetBodyTrackMask(int body_id, double body_params[],double track_params[],Body* body, Track* track){
    if (body_id == 0 ){
        Vacuum* body = new Vacuum;
        Vacuum::Track* track = new Vacuum::Track(track_params[0],track_params[1]);
    } else if(body_id == 1){
        ConstantDensity* body = new ConstantDensity(body_params[0],body_params[1]);
        ConstantDensity::Track* track = new ConstantDensity::Track(track_params[0],track_params[1]);
    } else if(body_id == 2){
        VariableDensity* body = new VariableDensity(body_params[0],body_params[1],body_params[2]);
        VariableDensity::Track* track = new VariableDensity::Track(track_params[0],track_params[1]);
    } else if(body_id == 3){
        Star* body = new Star(body_params[0],body_params[1],body_params[2]);
        Star::Track* track = new Star::Track(track_params[0],track_params[1]);
    } else if(body_id == 4){
        Earth* body = new Earth;
        Earth::Track* track = new Earth::Track(track_params[0],track_params[1],track_params[2]);
    } else if(body_id == 5){
        AGN* body = new AGN(body_params[0],body_params[1],body_params[2]);
        AGN::Track* track = new AGN::Track(track_params[0],track_params[1]);
    } else if(body_id == 6){
        Sun* body = new Sun;
        Sun::Track* track = new Sun::Track(track_params[0],track_params[1]);
    } else if(body_id == 7){
        SunASnu* body = new SunASnu;
        SunASnu::Track* track = new SunASnu::Track(track_params[0],track_params[1]);
    } else if(body_id == 8){
        EarthAtm* body = new EarthAtm;
        EarthAtm::Track* track = new EarthAtm::Track(track_params[0]);
    } else {
        printf("MaskError: Wrong body_id : %d",body_id);
    }
    
    return GSL_SUCCESS;
}

*/

BodyTrack* GetBodyTrackMask(int body_id, double body_params[],double track_params[]){
    BodyTrack* bodytrack = new BodyTrack;
    
    if (body_id == 0 ){
        Vacuum* body = new Vacuum;
        Vacuum::Track* track = new Vacuum::Track(track_params[0],track_params[1]);
        bodytrack->body = body;
        bodytrack->track = track;
    } else if(body_id == 1){
        ConstantDensity* body = new ConstantDensity(body_params[0],body_params[1]);
        ConstantDensity::Track* track = new ConstantDensity::Track(track_params[0],track_params[1]);
        bodytrack->body = body;
        bodytrack->track = track;
    } else if(body_id == 2){
        VariableDensity* body = new VariableDensity(body_params[0],body_params[1],body_params[2]);
        VariableDensity::Track* track = new VariableDensity::Track(track_params[0],track_params[1]);
        bodytrack->body = body;
        bodytrack->track = track;
    } else if(body_id == 3){
        Star* body = new Star(body_params[0],body_params[1],body_params[2]);
        Star::Track* track = new Star::Track(track_params[0],track_params[1]);
        bodytrack->body = body;
        bodytrack->track = track;
    } else if(body_id == 4){
        Earth* body = new Earth;
        Earth::Track* track = new Earth::Track(track_params[0],track_params[1],track_params[2]);
        bodytrack->body = body;
        bodytrack->track = track;
    } else if(body_id == 5){
        AGN* body = new AGN(body_params[0],body_params[1],body_params[2]);
        AGN::Track* track = new AGN::Track(track_params[0],track_params[1]);
        bodytrack->body = body;
        bodytrack->track = track;
    } else if(body_id == 6){
        Sun* body = new Sun;
        Sun::Track* track = new Sun::Track(track_params[0],track_params[1]);
        bodytrack->body = body;
        bodytrack->track = track;
    } else if(body_id == 7){
        SunASnu* body = new SunASnu;
        SunASnu::Track* track = new SunASnu::Track(track_params[0],track_params[1]);
        bodytrack->body = body;
        bodytrack->track = track;
    } else if(body_id == 8){
        EarthAtm* body = new EarthAtm;
        EarthAtm::Track* track = new EarthAtm::Track(track_params[0]);
        bodytrack->body = body;
        bodytrack->track = track;
    } else {
        printf("MaskError: Wrong body_id : %d",body_id);
    }
    
    return bodytrack;
}

/*/
 * NEUOSC MASK
 *
/*/

vector<double> CalNeuOscGSLMask(int ineu,double E,int body_id,double body_params[],double track_params[],int neutype,double paramlist[], double abs_error,double rel_error,int return_state,int optimization){
    // define stuff
    // configure oscillation parameters
    PhysConst* param = new PhysConst;
    SetOscParamMask(paramlist,param);
    // setting neutrino-antineutrino mode
    param->neutype = neutype;
    param->Refresh();    
    // calculate mixing matrix
    gsl_matrix_complex* fM2 = flavorM2(param);
    
    // configure body
    BodyTrack* bodytrack = GetBodyTrackMask(body_id,body_params,track_params);
    
    // oscillation probability vector definition
    int vec_size;
    if (return_state) {
        vec_size = 2*param->numneu;
    } else {
        vec_size = param->numneu;
    }
    double osc_prob[vec_size];
    
    #ifdef Mask_DEBUG
        cout << "Body Id : " << body_id << endl;
    #endif    
    
    CalNeuOscGSL(osc_prob,ineu,E,bodytrack->track,bodytrack->body,fM2,param,abs_error,rel_error,return_state,optimization);
    
    vector<double> prob;
    copy(&osc_prob[0],&osc_prob[vec_size],back_inserter(prob));
    
    return prob;
}


/*/
 * NEURHO MASK
 *
/*/

vector< vector<double> > CalNeuOscIntMask(int flux_type, double flux_params[],
                                        double Emin,double Emax, int Ediv,
                                        int body_id,double body_params[],double track_params[],
                                        int neutype,double paramlist[],
                                        double abs_error,double rel_error,
                                        int oscillation,int attenuation,int nc_inter,int cc_inter,
                                        int tau_regeneration,int muon,int electron,
                                        int basis){
    
    #ifdef CalNeuOscIntMask_DEBUG
        cout << "Setting up calculation" << endl;
    #endif
    
    // configure oscillation parameters
    PhysConst* param = new PhysConst;
    SetOscParamMask(paramlist,param);
    // setting neutrino-antineutrino mode
    param->neutype = neutype;
    // calculate mixing cosines/sines
    param->Refresh();   
    // configure body
    BodyTrack* bodytrack = GetBodyTrackMask(body_id,body_params,track_params);
        
    OscSetup* setup = new OscSetup;
    // setup the calculation
    setup->param = param;
    setup->E_range = logspace(Emin,Emax,Ediv);
    // setup interactions
    setup->oscillation = (oscillation != 0);
    setup->attenuation = (attenuation != 0);
    setup->nc_inter = (nc_inter != 0);
    setup->cc_inter = (cc_inter != 0);
    setup->tau_regeneration = (tau_regeneration != 0);
    setup->muon = (muon != 0);
    setup->electron = (electron != 0);
    // setup body
    setup->body = bodytrack->body;
    setup->track = bodytrack->track;
    // setup basis
    setup->basis = basis;

    VectorFluxState* initial_state;
    VectorFluxState* initial_state_aux = NULL;
    if (flux_type == 0){
        initial_state = InitSimpleState((int) flux_params[0],setup);
    } else if (flux_type == 1){
        initial_state = InitPowerLawFluxState((int) flux_params[0],flux_params[1],setup);
    } else if (flux_type == 2){
        initial_state = InitBodyFluxState(setup);
    } else if (flux_type == 9){
        int size = flux_params[0];
        
        FluxBox* Flux = new FluxBox(size);
        #ifdef CalNeuOscIntMask_DEBUG
            cout << "Begin Loading Generic Flux" << endl;
            cout << "size " << size << endl;
        #endif

        for(int i = 0; i < size; i++){
            Flux->E_values[i] = flux_params[1+size*0+i];
            
            Flux->nue_values[i] = flux_params[1+size*1+i];
            Flux->numu_values[i] = flux_params[1+size*2+i];
            Flux->nutau_values[i] = flux_params[1+size*3+i];
                        
            Flux->anue_values[i] = flux_params[1+size*4+i];
            Flux->anumu_values[i] = flux_params[1+size*5+i];
            Flux->anutau_values[i] = flux_params[1+size*6+i];
        }
        #ifdef CalNeuOscIntMask_DEBUG
            cout << "End Loading Generic Flux" << endl;        
        #endif     
        if( neutype != 2 ){
            initial_state = InitGeneralFlux(Flux,setup,neutype);
        } else {
            initial_state = InitGeneralFlux(Flux,setup,0);
            initial_state_aux = InitGeneralFlux(Flux,setup,1);    
        }     
    }
    
    VectorFluxState *initial_state_vector[2];
    initial_state_vector[0] = initial_state;
    initial_state_vector[1] = initial_state_aux;
    
    #ifdef CalNeuOscIntMask_DEBUG
        printf("Oscillation paramers :\n");
        printf("Emin : %g [eV]\n", Emin);
        printf("Emax : %g [eV]\n", Emax);
        cout << "Body Id : " << body_id << endl;
        printf("xini : %g [km]\n", setup->track->xini/param->km);
        printf("xend : %g [km]\n", setup->track->xend/param->km);
        printf("Square mass difference matrix \n");
        gsl_matrix_complex_fprintf(stdout,setup->flavorM2,"%g");
        
        printf("Get body information \n");
        cout << "name : " << setup->body->name << endl;
        printf("density : %g [gr/cm^3]\n",setup->body->density(setup->track));
    #endif

    #ifdef CalNeuOscIntMask_DEBUG
        cout << "Begin calculation" << endl;
    #endif

    VectorFluxState **final_state_vector = CalNeuOscInt(initial_state_vector,setup,abs_error,rel_error,true);
    
    #ifdef CalNeuOscIntMask_DEBUG
        cout << "End calculation" << endl;
    #endif    

    vector< vector<double> > state;
    vector<double> e_state;
    
    int e_size = setup->E_range.size();
    
    VectorFluxState* final_state = final_state_vector[0];
    VectorFluxState* final_state_aux = final_state_vector[1];
    
    for(int e = 0; e < e_size; e++){
        e_state.clear();
        e_state.push_back(setup->E_range[e]);
        for(int flavor = 0; flavor < setup->param->numneu; flavor++){
            e_state.push_back(
                TraceFlavorProject(flavor,final_state)->data[e]
                             );
        }
        state.push_back(e_state);
    }
    

    
    if (neutype == 2){
        for(int e = 0; e < e_size; e++){
            e_state.clear();
            e_state.push_back(setup->E_range[e]);
            for(int flavor = 0; flavor < setup->param->numneu; flavor++){
                e_state.push_back(
                    TraceFlavorProject(flavor,final_state_aux)->data[e]
                                 );
            }
            state.push_back(e_state);
        }
    }
    
    return state;
}


/*/
 * NEUSUN MASK
 *
/*/

vector< vector<double> > CalNeuOscSUNMask(int flux_type, double flux_params[],
                                        double Emin,double Emax, int Ediv,
                                        int body_id,double body_params[],double track_params[],
                                        int neutype,double paramlist[],
                                        double abs_error,double rel_error,
                                        int oscillation,int attenuation,int nc_inter,int cc_inter,
                                        int tau_regeneration,int muon,int electron,
                                        int basis){
    
    #ifdef CalNeuOscSUNMask_DEBUG
        cout << "Setting up calculation" << endl;
    #endif
    
    // configure oscillation parameters
    PhysConst* param = new PhysConst;
    SetOscParamMask(paramlist,param);
    // setting neutrino-antineutrino mode
    param->neutype = neutype;
    // calculate mixing cosines/sines
    param->Refresh();   
    // configure body
    BodyTrack* bodytrack = GetBodyTrackMask(body_id,body_params,track_params);
        
    OscSetup* setup = new OscSetup;
    // setup the calculation
    setup->param = param;
    setup->E_range = logspace(Emin,Emax,Ediv);
    // setup interactions
    setup->oscillation = (oscillation != 0);
    setup->attenuation = (attenuation != 0);
    setup->nc_inter = (nc_inter != 0);
    setup->cc_inter = (cc_inter != 0);
    setup->tau_regeneration = (tau_regeneration != 0);
    setup->muon = (muon != 0);
    setup->electron = (electron != 0);
    // setup body
    setup->body = bodytrack->body;
    setup->track = bodytrack->track;
    // setup basis
    setup->basis = basis;

    VectorFluxState* initial_state;
    VectorFluxState* initial_state_aux = NULL;
    if (flux_type == 0){
        initial_state = InitSimpleState((int) flux_params[0],setup);
    } else if (flux_type == 1){
        initial_state = InitPowerLawFluxState((int) flux_params[0],flux_params[1],setup);
    } else if (flux_type == 2){
        initial_state = InitBodyFluxState(setup);
    } else if (flux_type == 9){
        int size = flux_params[0];
        
        FluxBox* Flux = new FluxBox(size);
        #ifdef CalNeuOscSUNMask_DEBUG
            cout << "Begin Loading Generic Flux" << endl;
            cout << "size " << size << endl;
        #endif

        for(int i = 0; i < size; i++){
            Flux->E_values[i] = flux_params[1+size*0+i];
            
            Flux->nue_values[i] = flux_params[1+size*1+i];
            Flux->numu_values[i] = flux_params[1+size*2+i];
            Flux->nutau_values[i] = flux_params[1+size*3+i];
                        
            Flux->anue_values[i] = flux_params[1+size*4+i];
            Flux->anumu_values[i] = flux_params[1+size*5+i];
            Flux->anutau_values[i] = flux_params[1+size*6+i];
        }
        #ifdef CalNeuOscSUNMask_DEBUG
            cout << "End Loading Generic Flux" << endl;        
        #endif     
        if( neutype != 2 ){
            initial_state = InitGeneralFlux(Flux,setup,neutype);
        } else {
            initial_state = InitGeneralFlux(Flux,setup,0);
            initial_state_aux = InitGeneralFlux(Flux,setup,1);    
        }     
    }
    
    VectorFluxState *initial_state_vector[2];
    initial_state_vector[0] = initial_state;
    initial_state_vector[1] = initial_state_aux;
    
    #ifdef CalNeuOscSUNMask_DEBUG
        printf("Oscillation paramers :\n");
        printf("Emin : %g [eV]\n", Emin);
        printf("Emax : %g [eV]\n", Emax);
        cout << "Body Id : " << body_id << endl;
        printf("xini : %g [km]\n", setup->track->xini/param->km);
        printf("xend : %g [km]\n", setup->track->xend/param->km);
        
        printf("Get body information \n");
        cout << "name : " << setup->body->name << endl;
        printf("density : %g [gr/cm^3]\n",setup->body->density(setup->track));
    #endif

    #ifdef CalNeuOscSUNMask_DEBUG
        cout << "Begin calculation" << endl;
    #endif

    //VectorFluxState* final_state_vector;
    //final_state_vector = CalNeuOscSUN(initial_state_vector,setup,abs_error,rel_error);
    CalNeuOscSUN(initial_state_vector,setup,abs_error,rel_error);
    VectorFluxState** final_state_vector = initial_state_vector;
    
    #ifdef CalNeuOscSUNMask_DEBUG
        cout << "End calculation" << endl;
    #endif    

    vector< vector<double> > state;
    vector<double> e_state;
    
    int e_size = setup->E_range.size();
    
    VectorFluxState* final_state = final_state_vector[0];
    VectorFluxState* final_state_aux = final_state_vector[1];
    
    #ifdef CalNeuOscSUNMask_DEBUG
        cout << "Create vector states" << endl;
    #endif   
    
    for(int e = 0; e < e_size; e++){
        e_state.clear();
        e_state.push_back(setup->E_range[e]);
        cout << setup->param->numneu << endl;
        for(int flavor = 0; flavor < setup->param->numneu; flavor++){
            cout << TraceFlavorProject(flavor,final_state)->data[e] << endl;
            e_state.push_back(
                TraceFlavorProject(flavor,final_state)->data[e]
                             );
        }
        state.push_back(e_state);
    }
    
    if (neutype == 2){
        for(int e = 0; e < e_size; e++){
            e_state.clear();
            e_state.push_back(setup->E_range[e]);
            for(int flavor = 0; flavor < setup->param->numneu; flavor++){
                e_state.push_back(
                    TraceFlavorProject(flavor,final_state_aux)->data[e]
                                 );
            }
            state.push_back(e_state);
        }
    }
    
    #ifdef CalNeuOscSUNMask_DEBUG
        cout << "Send vector states" << endl;
    #endif   
    
    return state;
}
