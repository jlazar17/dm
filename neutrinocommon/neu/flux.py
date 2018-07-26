""" 
Author  : C.A. Arguelles
Date    : 25/MAR/2013

This package contains the flux definitions to be used
with the shareobject c++ interface.
"""
import neutrinocommon.neu.neuosc as no
import neutrinocommon.tools.generaltools as gt
import matplotlib.pyplot as plt
import numpy as np

class SimpleFlux():
    def __init__(self,ineu):
        self.id = 0
        self.fluxparams = [ineu]
        
class PowerLawFlux():
    def __init__(self,ineu,power):
        self.id = 1
        self.fluxparams = [ineu,power]
        
class GeneralFlux():
    def __init__(self, neu_eq_aneu = False):
        self.id = 9
        self.neu_eq_aneu = neu_eq_aneu
        self.E_range = []
        
        self.nue_values = []
        self.numu_values = []
        self.nutau_values = []
        
        self.anue_values = []
        self.anumu_values = []
        self.anutau_values = []
        
    def load(self):
        self.size = len(self.E_range)
        
        if self.neu_eq_aneu :
            self.anue_values = self.nue_values
            self.anumu_values = self.numu_values
            self.anutau_values = self.nutau_values
            
        if len(self.nue_values) == 0 :
            self.nue_values = [0.0]*self.size

        if len(self.numu_values) == 0 :
            self.numu_values = [0.0]*self.size

        if len(self.nutau_values) == 0 :
            self.nutau_values = [0.0]*self.size
            
        if len(self.anue_values) == 0 :
            self.anue_values = [0.0]*self.size

        if len(self.anumu_values) == 0 :
            self.anumu_values = [0.0]*self.size

        if len(self.anutau_values) == 0 :
            self.anutau_values = [0.0]*self.size
        
        self.fluxparams = [self.size] + self.E_range + self.nue_values + self.numu_values + self.nutau_values + \
                            self.anue_values + self.anumu_values + self.anutau_values
        
        
def FluxToFile(flux,filepath):
    data = []
    
    for i,e in enumerate(flux.E_range):
        data.append([e,
                     flux.nue_values[i],flux.numu_values[i],flux.nutau_values[i],
                     flux.anue_values[i],flux.anumu_values[i],flux.anutau_values[i]])
        
    file = open(filepath,'w')
    gt.quickprint(file,data)
    file.close()


def FileToFlux(flux,filepath):
    file = open(filepath,'r')
    data = gt.quickread(file)
    file.close()
    
    flux.E_range = [ dat[0] for dat in data]
    
    flux.nue_values = [ dat[1] for dat in data]
    flux.numu_values = [ dat[2] for dat in data]
    flux.nutau_values = [ dat[3] for dat in data]
    
    flux.anue_values = [ dat[4] for dat in data]
    flux.anumu_values = [ dat[5] for dat in data]
    flux.anutau_values = [ dat[6] for dat in data]
    
def FileToFMFlux(flux,filepath, basis = 'flavor'):
    file = open(filepath,'r')
    data = gt.quickread(file)
    file.close()
    
    flux.E_range = [ dat[0] for dat in data]
    
    numneu = (len(data[0])-1)/4
    
    if basis == 'flavor':
        flux.nue_values = [ dat[1] for dat in data]
        flux.numu_values = [ dat[2] for dat in data]
        flux.nutau_values = [ dat[3] for dat in data]
        
        flux.anue_values = [ dat[4] for dat in data]
        flux.anumu_values = [ dat[5] for dat in data]
        flux.anutau_values = [ dat[6] for dat in data]
    elif basis == 'mass':
        flux.nue_values = [ dat[1+2*numneu] for dat in data]
        flux.numu_values = [ dat[2+2*numneu] for dat in data]
        flux.nutau_values = [ dat[3+2*numneu] for dat in data]
        
        flux.anue_values = [ dat[4+2*numneu] for dat in data]
        flux.anumu_values = [ dat[5+2*numneu] for dat in data]
        flux.anutau_values = [ dat[6+2*numneu] for dat in data]
    else :
        quit()
    
def FileToRhoFlux(filepath):
    
    datafile = open(filepath,'r')
    data_array = gt.quickread(datafile)
    datafile.close()
    
    E_range = []
    rho = []
    arho = []
    for data in data_array:
        E_range.append(data[0])
        
        numarray = (len(data)-1)/2 + 1
        
        data_neu = data[1:numarray]
        data_neu_complex = [ complex(data_neu[i],data_neu[i+1]) for i in range(0,len(data_neu)-1,2) ] 
        data_aneu = data[numarray:len(data)]
        data_aneu_complex = [ complex(data_aneu[i],data_aneu[i+1]) for i in range(0,len(data_aneu)-1,2) ] 
        
        numneu = int(np.sqrt(len(data_neu_complex)))

        neu = np.zeros([numneu,numneu],complex)
        aneu = np.zeros([numneu,numneu],complex)
    
        for i in range(numneu):
            for j in range(numneu):
                neu[i,j] = data_neu_complex[numneu*i+j]
                aneu[i,j] = data_aneu_complex[numneu*i+j]
                
        rho.append(neu)
        arho.append(aneu)
        
    return E_range,rho,arho

def RhoFluxToFlux(E_range,rho,arho,
                  param = None,
                  complete_decoherence = False):
    
    if complete_decoherence and param == None:
        print "ERROR: Need to submit PhysicsConstants in order to process."
        quit()
        
    numneu = int(np.sqrt(rho[0].size))
    
    # construct projector matrices
    pi_array = []
    for i in range(numneu):
        # build flavor projector
        pi = np.zeros([3,3],complex)
        pi[i,i] = complex(1.0,0.0)
        if complete_decoherence :
            # convert to mass projector
            PMNS = no.MixMatrix(param)
            pi = np.dot(PMNS.U,np.dot(pi,PMNS.UCT))
        
        pi_array.append(pi)
    
    flux = GeneralFlux()
    
    flux.E_range = E_range
    
    for neu in rho:
        flux.nue_values.append(np.trace(np.dot(pi_array[0],neu)).real)
        flux.numu_values.append(np.trace(np.dot(pi_array[1],neu)).real)
        flux.nutau_values.append(np.trace(np.dot(pi_array[2],neu)).real)
        
    for aneu in arho:
        flux.anue_values.append(np.trace(np.dot(pi_array[0],aneu)).real)
        flux.anumu_values.append(np.trace(np.dot(pi_array[1],aneu)).real)
        flux.anutau_values.append(np.trace(np.dot(pi_array[2],aneu)).real)
        
    flux.load()
    
    return flux

### PLOT TOOLS ###

def PlotFlux(flux,e_power = 0,scale = 'loglog',
                xlim = None,
                ylim = None):
    fig = plt.figure(figsize = (14,7))
    plt.rc('font', family='serif', size=20)
    plt.rc('legend', fontsize = 15)
    ax1 = plt.subplot(121)

    plt.xlabel(r'$E_\nu [GeV]$')
    if e_power != 0 :
        plt.ylabel(r'$F_\nu E^'+str(e_power)+' $')
    else : 
        plt.ylabel(r'$F_\nu $')
        
    if scale == 'loglog':
        plt.loglog()
    elif scale == 'linlog':
        plt.semilogy()
    elif scale == 'loglin':
        plt.semilogx()
        
    if xlim != None :
        plt.xlim(xlim)
    if ylim != None :
        plt.ylim(ylim)
        
    ax2 = plt.subplot(122)
    
    plt.xlabel(r'$E_\nu [GeV]$')
    if e_power != 0 :
        plt.ylabel(r'$F_\nu E^'+str(e_power)+' $')
    else : 
        plt.ylabel(r'$F_\nu $')

    if scale == 'loglog':
        plt.loglog()
    elif scale == 'linlog':
        plt.semilogy()
    elif scale == 'loglin':
        plt.semilogx()
        
    if xlim != None :
        plt.xlim(xlim)
    if ylim != None :
        plt.ylim(ylim)
        
    E_range_GeV = [e*1.0e-9 for e in flux.E_range]

    nue   = [f*flux.E_range[i]**e_power for i,f in enumerate(flux.nue_values)]
    numu  = [f*flux.E_range[i]**e_power for i,f in enumerate(flux.numu_values)]
    nutau = [f*flux.E_range[i]**e_power for i,f in enumerate(flux.nutau_values)]
    
    anue   = [f*flux.E_range[i]**e_power for i,f in enumerate(flux.anue_values)]
    anumu  = [f*flux.E_range[i]**e_power for i,f in enumerate(flux.anumu_values)]
    anutau = [f*flux.E_range[i]**e_power for i,f in enumerate(flux.anutau_values)]
    
    ax1.plot(E_range_GeV,nue,label = r'$ \nu_e $',lw = 3)
    ax1.plot(E_range_GeV,numu,label = r'$ \nu_\mu $',lw = 3)
    ax1.plot(E_range_GeV,nutau,label = r'$ \nu_\tau $',lw = 3)
    ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fancybox = True)
    
    ax2.plot(E_range_GeV,anue,label = r'$ \bar{\nu}_e $',lw = 3)
    ax2.plot(E_range_GeV,anumu,label = r'$ \bar{\nu}_\mu $',lw = 3)
    ax2.plot(E_range_GeV,anutau,label = r'$ \bar{\nu}_\tau $',lw = 3)
    ax2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fancybox = True)
    
    plt.subplots_adjust(wspace = 0.5)
    
def PlotCombineFlux(flux,e_power = 0,scale = 'loglog',
                    xlim = None, xscale = 1.0,
                    ylim = None, yscale = 1.0,
                    ylabel = None, xlabel = None,
                    textloc = None, text = None,
                    lloc = None, savepath = None,
                    add = False):
    
    fig = plt.figure(figsize = (7,7))
    plt.rc('font', family='serif', size=20)
    plt.rc('legend', fontsize = 15)

    if xlabel != None:
        plt.xlabel(xlabel)
    else:
        plt.xlabel(r'$E_\nu [GeV]$')
    
    if ylabel != None:
        plt.ylabel(ylabel)
    else :
        if e_power != 0 :
            plt.ylabel(r'$F_\nu E^'+str(e_power)+' $')
        else : 
            plt.ylabel(r'$F_\nu $') 
        
    if scale == 'loglog':
        plt.loglog()
    elif scale == 'linlog':
        plt.semilogy()
    elif scale == 'loglin':
        plt.semilogx()
        
    if xlim != None :
        plt.xlim(xlim)
    if ylim != None :
        plt.ylim(ylim)
        
    E_range_GeV = [e*1.0e-9*xscale for e in flux.E_range]
    
    nue   = [f*yscale*flux.E_range[i]**e_power for i,f in enumerate(flux.nue_values)]
    numu  = [f*yscale*flux.E_range[i]**e_power for i,f in enumerate(flux.numu_values)]
    nutau = [f*yscale*flux.E_range[i]**e_power for i,f in enumerate(flux.nutau_values)]
    
    anue   = [f*yscale*flux.E_range[i]**e_power for i,f in enumerate(flux.anue_values)]
    anumu  = [f*yscale*flux.E_range[i]**e_power for i,f in enumerate(flux.anumu_values)]
    anutau = [f*yscale*flux.E_range[i]**e_power for i,f in enumerate(flux.anutau_values)]
    
    if add :
        num_p = len(nue)
        nue_add   = [ nue[i] + anue[i] for i in range(num_p) ]
        numu_add  = [ numu[i] + anumu[i] for i in range(num_p) ]
        nutau_add = [ nutau[i] + anutau[i] for i in range(num_p) ]
    
    cols=['#29A2C6','#FF6D31','#FFCB18','#73B66B']
    
    if not add :
        plt.plot(E_range_GeV,nue,label = r'$ \nu_e $',lw = 3, linestyle = 'solid', color = cols[0])
        plt.plot(E_range_GeV,numu,label = r'$ \nu_\mu $',lw = 3, linestyle = 'solid', color = cols[1])
        plt.plot(E_range_GeV,nutau,label = r'$ \nu_\tau $',lw = 3, linestyle = 'solid', color = cols[3])
        
        plt.plot(E_range_GeV,anue,label = r'$ \bar{\nu}_e $',lw = 3, linestyle = 'dashed', color = cols[0])
        plt.plot(E_range_GeV,anumu,label = r'$ \bar{\nu}_\mu $',lw = 3, linestyle = 'dashed', color = cols[1])
        plt.plot(E_range_GeV,anutau,label = r'$ \bar{\nu}_\tau $',lw = 3, linestyle = 'dashed', color = cols[3])
    else :
        plt.plot(E_range_GeV,nue_add,label = r'$ \nu_e + \bar{\nu}_e $',lw = 3, linestyle = 'solid', color = cols[0])
        plt.plot(E_range_GeV,numu_add,label = r'$ \nu_\mu + \bar{\nu}_\mu$',lw = 3, linestyle = 'solid', color = cols[1])
        plt.plot(E_range_GeV,nutau_add,label = r'$ \nu_\tau + \bar{\nu}_\tau$',lw = 3, linestyle = 'solid', color = cols[3])    
    
    if text != None and textloc != None :
        plt.annotate(text,xy = textloc, xytext = textloc)
    
    if lloc != None :
        plt.legend(loc = lloc, fancybox = True)
    else :
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fancybox = True)

    
    plt.gcf().subplots_adjust(left=0.17)

    if savepath != None:
        plt.savefig(savepath, dpi = 600)