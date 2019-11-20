#!/usr/bin/env python3

#Example to describe the main steps required to define the model and inpute paramaters
#and solve the boltzmann equations

#First tell the system where to find the modules:
import os
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)




def main(parameterFile,outputFile,showPlot=True):
    """
    
    Main code to define the BSM contents and properties and solve the Boltzmann equations
    
    :param parameterFile: Path to the file defining the main model parameters
    :param outputFile: Path to the output file. If None, no results will be written.
    :param showPlot: If True, will show a simple plot for the evolution of the energy densities
    
    """

    from pyCode.component import Component
    from pyCode.boltzSolver import BoltzSolution
    from pyCode.AuxDecays import DecayList, Decay
    import numpy as np
    from scipy.interpolate import interp1d
    
    lamb = 2.6e-7
#     lamb = 0.17
    decays = DecayList()
    decayToDM = Decay(instate='Mediator',fstates=['DM','radiation'],br=1.)
    decays.addDecay(decayToDM)
    decays.Xfraction = 0.1
    decays.width = 2.49e-15*(lamb/4.3e-7)**2
    

    #Get the model parameters (or define them here):
    TRH = 1e4
    TF = 1e-3
    
    def nEQbottom(T):
        Zeta3 = 1.20206
        return 3*2*(3./4.)*Zeta3*T**3/np.pi**2

    def nEQgluon(T):
        Zeta3 = 1.20206
        return 8*2*Zeta3*T**3/np.pi**2
    
    #Annihilation rate for mediator
    data = np.genfromtxt('./width_and_medxs.dat',skip_header=5)
    dataR = np.genfromtxt('./sigmav_conversion_bchi-sbotg_500_510_1.dat',skip_header=6,usecols=(0,4))
    conv = 0.8579e17

    sLog = lambda x: interp1d(data[:,0],np.log(data[:,1]*conv),
                        fill_value='extrapolate',bounds_error=False)(x)
    cRateLog = lambda x: interp1d(dataR[:,0],np.log(dataR[:,1]*conv*(lamb)**2),
                        fill_value='extrapolate',bounds_error=False)(x)

    #Conversion rates for DM and mediator: 
    dofDM = -2 #Number of DM degrees of freedom (Majorana fermion)
    dofMed = 6 #Number of Mediator degrees of freedom (complex colored scalar)

    @np.vectorize
    def sigmaVJan(T):
        x = 500./T
        if x > data[:,0].max():
            return 0.
        sF = sLog(x)
        return np.exp(sF)

    @np.vectorize
    def cRateDMJan(T):
        x = 500./T
        if x > dataR[:,0].max():
            return 0.
        sF = cRateLog(x)
        return 2*nEQbottom(T)*np.exp(sF)

    @np.vectorize
    def cRateMedJan(T):
        x = 510./T
        if x > dataR[:,0].max():
            return 0.
        sF = cRateLog(x)
        return 2*nEQgluon(T)*np.exp(sF)*abs(dofMed/dofDM)*np.exp((510.-500.)/T)

    
    #Define the components to be evolved and their properties:    
    dm = Component(label='DM',Type='thermal',dof=dofDM,
                   mass=500.
                    ,sigVij=lambda T,other: 1e-10*sigmaVJan(T)
                    ,sigVii=lambda T: 1e-10*sigmaVJan(T)
                    ,cRate=lambda T,other: cRateDMJan(T)
                   )
    mediator = Component(label='Mediator',Type='thermal',dof=dofMed,
                   mass=510.,decays=decays,
                   sigVii=sigmaVJan
                    ,cRate=lambda T,other: cRateMedJan(T)
                   )
    compList = [dm,mediator]
    
    #Evolve the equations from TR to TF
    solution = BoltzSolution(compList,TRH)
    solved = solution.EvolveTo(TF,npoints=5000)
    if not solved:
        return
    
    #Print summary
    if outputFile:
        if os.path.isfile(outputFile):
            os.remove(outputFile)
        solution.printSummary(outputFile)
        solution.printData(outputFile)
    solution.printSummary()
    

if __name__ == "__main__":

    import argparse    
    ap = argparse.ArgumentParser( description=
            "Evolve Boltzmann equations for a simple non-thermal DM scenario" )
    ap.add_argument('-p', '--parameterFile', 
            help='name of parameter file, where most options are defined', default = 'parameters.ini')
    ap.add_argument('-o', '--outputFile', 
            help='name of output file (optional argument). If not define, no output will be saved', default=None)
    ap.add_argument('-P', '--plotResult', help='show simple plot for the evolution of densities',
            action='store_true',default = False)    
    
    args = ap.parse_args()
    main(args.parameterFile,args.outputFile,args.plotResult)
