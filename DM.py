import numpy as np
import scipy as sp
import scipy.special as spe
import scipy.interpolate as interpolate

def A(beta,vel,velChar):
    return np.sqrt(beta)*vel/velChar

def AUpPlus(A,a,etaUp):
    return A*np.sqrt(1+a) + etaUpHat

def AUpMinus(A,a,etaUp):
    return A*np.sqrt(1+a) - etaUpHat

def ADnPlus(A,b,etaDn):
    return A*np.sqrt(1+b) + etaDnHat

def ADnMinus(A,b,etaDn):
    return A*np.sqrt(1+b) - etaDnHat

# Calculates term outside parenthesis in eq 2.18
def gouldCoefficient(sigma0,dmRho,massSun,velChar,Q,ep,eta,a):
    return sigma0*dmRho*massSun*velChar*Q**2*ep/(2*eta*a)

# Calculates first term in parenthesis in eq 2.18
def gouldFirstTerm(ACenter,ASurface,a,eta):
    tmp = 2*np.exp(-a*eta**2)spe.erf(eta)/(1+a)
    result = tmp - np.exp(-a*eta)/((ACenter**2-ASurface**2)*np.power((1+a),3./2.))
    return result

# Calculates second term in parenthesis in eq 2.18
def gouldSecondTerm(A,a,b,etaUp):
    APlus  = AUpPlus(A,a,etaUp)
    AMinus = AUpMinus(A,s,etaUp)
    tmp    = (APLus*AMinus-1./2-(1+a)/(a-b))(spe.erf(APlus)-spe.erf(AMinus))
    result = tmp + (AMinus*np.exp(-APlus**2)-APlus*np.exp(-AMinus**2))/np.pi
    return result

# Calculates third term in parenthesis in eq 2.18
def gouldThirdTerm(A,ASurface,ACenter,a,b,etaDn):
    APlus  = ADnPlus(A,a,etaDn)
    AMinus = ADnMinus(A,s,etaDn)
    tmp = (2*spe.erf(eta)-spe.erf(APlus)+spe.erf(AMinus))*np.exp((b-a)A**2)
    result = tmp*np.exp(-b*eta**2)/((a-b)(ACenter**2-ASurface**2)np.sqrt(1+b))
    return result

def DMSunCaptureRate(dmMass, dmCs, param):

    """ Calculates DM capture rate in the sun using Gould's formula.

    Ref: Cosmological density of WIMPs from solar and terresterial annihilations. A. Gould.
    The Astrophysical Journal, 388:338-344, 1992 April
    """

    elements        = {0:'H',1:'He',2:'C',3:'N',4:'O',5:'Ne',6 :'Mg',7:'Si',8:'S',9:'Fe'}              # elements in sun
    massNum         = [1.0,4.0,12.0,14.0,16.0,20.0,24.0,28.0,32.0,56.0]                                # A : mass number
    massGrPerMol    = [1.0079,4.0026,12.0107,14.0067,15.9994,20.1797,24.3051,28.0855,32.0655,55.8452]  # gr mol^-1
    eps             = [0.772,0.209,3.87e-3,9.4e-4,8.55e-3,1.51e-3,7.39e-4,8.13e-4,4.65e-4,1.46e-3]     # relative abundances in Sun from tbl 8 Jungman p.299
    n = len(elements)

    # input data
    mass_eV         = [ m*param.gr/param.Na for m in massGrPerMol ]
    atomicRadius    = [ 1.2*param.fermi*np.power(A, 1./3.)/param.Na for A in massNum ]
    energyElement   = [ 3. / (2*mass_eV[i]*atomicRadius[i]) for i in range(n) ]           # energy defined in Gould eq 2.5
    # is this still the best number to use
    dmRho           = 0.3*param.GeV/param.cm**3                                          # dm density Ref : arXiv 1105.6339
    massSun         = 1.9891e30*param.kg
    eta             = 1.0

    """ Gould's velocity definitions. 
    Are these still good numbers?
    """
    velCenter       = 1354*param.km/param.sec  # Gould eq 2.14
    velSurface      = 795*param.km/param.sec   # Gould eq 2.14
    velRot          = 220.*param.km/param.sec  # Gould page 341 below eq 19
    velChar         = velRot/eta

    # Jungman velocity definitions
    velBar          = 220.*param.km/param.sec
    velStar         = np.sqrt(3.0/2.0)*velBar/eta

    # special assumptions
    q = [ m for m in mass_eV ]  # arbitrary parameter from eq 2.4

    # define Gould's auxilary variables
    betaPlus        = [ 4.*dmMass*m/((dmMass+m)**2) for m in mass_eV ]
    betaMinus       = [ 4.*dmMass*m/((dmMass-m)**2) for m in mass_eV ]
    #sigma0          = [ dmCs for i in range(n) ]
    sigma0          = dmCs
    a               = [ dmMass*velChar**2/(2.*e) for e in energyElement ]
    b               = [ betaPlus[i]*a[i] for i in range(n) ]
    etaUp           = [ eta/np.sqrt(1+x) for x in a ]
    etaDn           = [ eta/np.sqrt(1+x) for x in b ]
    #A               = lambda beta, vel : np.sqrt(beta)*vel/velChar
    #AUpHatPlus      = lambda A, a : A*np.sqrt(1+a) + etaUpHat
    #AUpHatMinus     = lambda A, a : A*np.sqrt(1+a) - etaUpHat
    #ADnHatPlus      = lambda A, b : A*np.sqrt(1+b) + etaDnHat
    #ADnHatMinus     = lambda A, b : A*np.sqrt(1+b) - etaDnHat
    ACenter         = [ A(x,velCenter,velChar) for x in betaMinus ]
    ASurface        = [ A(x,velSurface,velChar) for x in betaMinus ]

    coefs           = [ gouldCoefficient(sigma0,dmRho,massSun,velChar,q[i],eps[i],eta,a[i]) \
                        for i in range(n) ]
    T1              = [ gouldFirstTerm(ACenter[i],ASurface[i],a[i],etaUp[i] \
                        for i in range(n) ]
    T2              = [ gouldSecondTerm(ACenter[i],a[i],b[i],etaUp[i] \
                        - gouldSecondTerm(ASurface[i],a[i],b[i],etaUp[i] \
                        for i in range(n) ]
    T3              = [ gouldThirdTerm(ACenter[i],ASurface[i],ACenter[i],a[i],b[i],etaDn[i]) \
                        - gouldThirdTerm(ASurface[i],ASurface[i],ACenter[i],a[i],b[i],etaDn[i]) \
                        for i in range(n) ]
    rate            = np.sum( [coefs[i]*T1[i]*(T2[i]+T3[i]) for i in range(n) ] )
    return rate
         
