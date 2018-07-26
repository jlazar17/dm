from libc.stdlib cimport malloc,free
from libcpp cimport bool
from libcpp.vector cimport vector

# creating link
cdef extern from "interface.h":
    vector[double] CalNeuOscGSLMask(int ineu,double E,
                                    int body_id,double body_params[],double track_params[],
                                    int neutype,double paramlist[],
                                    double abs_error,double rel_error,int return_state,int optimization)
    
    vector[vector[double]] CalNeuOscIntMask(int flux_type, double flux_params[],
                                        double Emin,double Emax, int div,
                                        int body_id,double body_params[],double track_params[],
                                        int neutype,double paramlist[],
                                        double abs_error,double rel_error,
                                        int oscillation,int attenuation,int nc_inter, int cc_inter,
                                        int tau_regeneration, int muon, int electron,
                                        int basis)
    vector[vector[double]] CalNeuOscSUNMask(int flux_type, double flux_params[],
                                        double Emin,double Emax, int div,
                                        int body_id,double body_params[],double track_params[],
                                        int neutype,double paramlist[],
                                        double abs_error,double rel_error,
                                        int oscillation,int attenuation,int nc_inter, int cc_inter,
                                        int tau_regeneration, int muon, int electron,
                                        int basis)

# aux interfaces
    
cdef double* MapPListToCList(list Plist):
    cdef int length = len(Plist)
    cdef int i
    
    cdef double *Clist = <double*>malloc(length*sizeof(double))
    for i in range(length):
        Clist[i] = Plist[i]
        
    return Clist

cdef list MapCListToPList(double* Clist,int length):
    cdef list Plist = []
    print "got here"
    print Clist[1]
    Plist.append(<double>Clist[1])
    print "died here"
    for i in range(length):
        Plist.append(<double>Clist[i])
    print "passed"
    return Plist

cdef double* OscParamTranslator(param):
    # get the oscillation parameters from the param object
    # number of flavors
    numneu = [param.numneu]
    # mixing angles
    mixangles = [param.th12,param.th13,param.th23,
                 param.th14,param.th24,param.th34,
                 param.th15,param.th25,param.th35,
                 param.th45,param.th16,param.th26,
                 param.th36,param.th46,param.th56]
    # square mass differences
    massdiff = [param.dm21sq,param.dm31sq,param.dm41sq,param.dm51sq,param.dm51sq]
    # cp-phases
    cpphases = [param.delta1,param.delta2,param.delta3]
    # neutrino species
    neutrinos = [param.electron,param.muon,param.tau,
                 param.sterile1,param.sterile2,param.sterile3]
    
    # join and return
    paramlist = numneu+mixangles+massdiff+cpphases+neutrinos
    return MapPListToCList(paramlist)

# main interface

cpdef CalNeuOscCPP(ineu,E,track,body,param,
                   NCCC_int = False,
                   return_state = False,complete_decoherence = False,
                   bool survival_prob = False,double E_threshold = 1.0,
                   bool step_osc_prob = False,bool only_osc_prob = True,
                   bool optimization = True,
                   double abs_error = 1.0e-8,double rel_error = 1.0e-8):
    """
    Call the CPP-RK implementations and solves Schrodinger equation.
    
    @type ineu : int
    @param ineu : initial flavor
    @type E : float
    @param E : neutrino energy
    @type Track : track
    @param Track : trayectory in the body
    @type Body  : body
    @param Body : body
    @type param : PhysConstant
    @param param : Oscillation parameters

    @rtype  : array
    @return : Oscillation Probaility
    
    """
    # body and track   
    cdef int body_id = body.id

    cdef double* body_params = MapPListToCList(body.bodyparams)
    cdef double* track_params = MapPListToCList(track.trackparams)
    
    # oscillation setup
    cdef double* paramlist = OscParamTranslator(param)
    
    cdef int neutype = 0
    if param.neutype == "antineutrino":
        neutype = 1
        
    cdef int return_state_int = 0
    if return_state :
        return_state_int = 1
    cdef int optimization_int = 0
    if optimization:
        optimization_int = 1

    # sent to CPP
    osc_prob = CalNeuOscGSLMask(ineu,E,body_id,body_params,track_params,
                     neutype,paramlist,
                     abs_error,rel_error,return_state_int,optimization_int)
    
    # free memory
    free(body_params)
    free(track_params)
    free(paramlist)
    # return

    return osc_prob


cpdef CalNeuIntCPP(flux, double Emin, double Emax, int div,
                   track,body,
                   param,
                   bool oscillation = False, bool attenuation = False,
                   bool nc_inter = False, bool cc_inter = False, bool tau_regeneration = False,
                   bool muon = False, bool electron = False,
                   double abs_error = 1.0e-20,double rel_error = 1.0e-20,
                   basis = 0):
    """
    Call the CPP-RK implementation of neutrino propagation in the
    density matrix formalism including oscillations,
    NC-interactions, CC-inteactions, and tau regeneration.
    
    @type flux : Flux
    @param flux : initial flux type
    @type Emin : float
    @param Emin : minimun neutrino energy
    @type Emax : float
    @param Emax : maximun neutrino energy
    @type div : int
    @param div : number of divisions
    @type Track : track
    @param Track : trayectory in the body
    @type Body  : body
    @param Body : body
    @type param : PhysConstant
    @param param : Oscillation parameters

    @rtype  : array
    @return : oscillation probability for several energies
    
    """
    # flux
    cdef int flux_type_id = flux.id
    cdef double* flux_params = MapPListToCList(flux.fluxparams)
    
    # body and track   
    cdef int body_id = body.id

    cdef double* body_params = MapPListToCList(body.bodyparams)
    cdef double* track_params = MapPListToCList(track.trackparams)
    
    # oscillation setup
    cdef double* paramlist = OscParamTranslator(param)
    
    cdef int neutype_int = 0
    if param.neutype == "antineutrino":
        neutype_int = 1
    elif param.neutype == "both":
        neutype_int = 2
    
    # options
    cdef int oscillation_int = 0
    if oscillation :
        oscillation_int = 1
    cdef int attenuation_int = 0
    if attenuation:
        attenuation_int = 1
    cdef int nc_inter_int = 0
    if nc_inter:
        nc_inter_int = 1
    cdef int cc_inter_int = 0
    if cc_inter:
        cc_inter_int = 1
    cdef int tau_regeneration_int = 0
    if tau_regeneration:
        tau_regeneration_int = 1  
    cdef int muon_int = 0
    if muon:
        muon_int = 1
    cdef int electron_int = 0
    if electron:
        electron_int = 1
        
    # basis
    cdef int basis_int = basis
    
    # sent to CPP
    e_osc_prob = CalNeuOscIntMask(flux_type_id,flux_params,
                                Emin,Emax,div,
                                body_id,body_params,track_params,
                                neutype_int,paramlist,
                                abs_error,rel_error,
                                oscillation_int,attenuation_int,
                                nc_inter_int,cc_inter_int,
                                tau_regeneration_int,muon_int,electron_int,
                                basis_int)
    
    # free memory
    free(body_params)
    free(track_params)
    free(flux_params)
    free(paramlist)
    # return

    if neutype_int != 2 : 
        return e_osc_prob
    else :
        e_osc_prob_neu = [ p for i,p in enumerate(e_osc_prob) if i <= div + 1]
        e_osc_prob_aneu = [ p for i,p in enumerate(e_osc_prob) if i > div + 1]
        return [e_osc_prob_neu,e_osc_prob_aneu]
    
cpdef CalNeuSUNCPP(flux, double Emin, double Emax, int div,
                   track,body,
                   param,
                   bool oscillation = False, bool attenuation = False,
                   bool nc_inter = False, bool cc_inter = False, bool tau_regeneration = False,
                   bool muon = False, bool electron = False,
                   double abs_error = 1.0e-20,double rel_error = 1.0e-20,
                   basis = 0):
    """
    Call the CPP-RK implementation of neutrino propagation in the
    density matrix formalism including oscillations,
    NC-interactions, CC-inteactions, and tau regeneration.
    
    @type flux : Flux
    @param flux : initial flux type
    @type Emin : float
    @param Emin : minimun neutrino energy
    @type Emax : float
    @param Emax : maximun neutrino energy
    @type div : int
    @param div : number of divisions
    @type Track : track
    @param Track : trayectory in the body
    @type Body  : body
    @param Body : body
    @type param : PhysConstant
    @param param : Oscillation parameters

    @rtype  : array
    @return : oscillation probability for several energies
    
    """
    # flux
    cdef int flux_type_id = flux.id
    cdef double* flux_params = MapPListToCList(flux.fluxparams)
    
    # body and track   
    cdef int body_id = body.id

    cdef double* body_params = MapPListToCList(body.bodyparams)
    cdef double* track_params = MapPListToCList(track.trackparams)
    
    # oscillation setup
    cdef double* paramlist = OscParamTranslator(param)
    
    cdef int neutype_int = 0
    if param.neutype == "antineutrino":
        neutype_int = 1
    elif param.neutype == "both":
        neutype_int = 2
    
    # options
    cdef int oscillation_int = 0
    if oscillation :
        oscillation_int = 1
    cdef int attenuation_int = 0
    if attenuation:
        attenuation_int = 1
    cdef int nc_inter_int = 0
    if nc_inter:
        nc_inter_int = 1
    cdef int cc_inter_int = 0
    if cc_inter:
        cc_inter_int = 1
    cdef int tau_regeneration_int = 0
    if tau_regeneration:
        tau_regeneration_int = 1  
    cdef int muon_int = 0
    if muon:
        muon_int = 1
    cdef int electron_int = 0
    if electron:
        electron_int = 1
        
    # basis
    cdef int basis_int = basis
    
    # sent to CPP
    e_osc_prob = CalNeuOscSUNMask(flux_type_id,flux_params,
                                Emin,Emax,div,
                                body_id,body_params,track_params,
                                neutype_int,paramlist,
                                abs_error,rel_error,
                                oscillation_int,attenuation_int,
                                nc_inter_int,cc_inter_int,
                                tau_regeneration_int,muon_int,electron_int,
                                basis_int)
    
    # free memory
    free(body_params)
    free(track_params)
    free(flux_params)
    free(paramlist)
    # return

    if neutype_int != 2 : 
        return e_osc_prob
    else :
        e_osc_prob_neu = [ p for i,p in enumerate(e_osc_prob) if i <= div + 1]
        e_osc_prob_aneu = [ p for i,p in enumerate(e_osc_prob) if i > div + 1]
        return [e_osc_prob_neu,e_osc_prob_aneu]

# tau decay extensions and links

cdef extern from "taudecay.h":
    double TauDecayToAll(double E_tau,double E_nu)
    double TauDecayToLepton(double E_tau,double E_nu)
cpdef double TauDecayToAllCPP(double E_tau,double E_nu):
    """
    Invokes the formulaes implemented in C++ for tau decay.
    
    Returns dn/dz where z = E_nu/E_tau
    
    @type E_tau : float
    @param E_tau : Tau energy [GeV]
    @type E_nu : float
    @param E_nu : Neutrino energy [GeV]

    @rtype  : float
    @return : dn/dz(E_tau,E_lep)
    
    """
    return TauDecayToAll(E_tau,E_nu)

cpdef double TauDecayToLeptonCPP(double E_tau,double E_nu):
    """
    Invokes the formulaes implemented in C++ for tau decay.
    
    Returns dn/dz where z = E_nu/E_tau
    
    @type E_tau : float
    @param E_tau : Tau energy [GeV]
    @type E_nu : float
    @param E_nu : Neutrino energy [GeV]

    @rtype  : float
    @return : dn/dz(E_tau,E_lep)
    
    """
    return TauDecayToLepton(E_tau,E_nu)