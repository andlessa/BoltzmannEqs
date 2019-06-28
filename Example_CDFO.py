#!/usr/bin/env python3

#Example to describe the main steps required to define the model and inpute paramaters
#and solve the boltzmann equations

#First tell the system where to find the modules:
import os
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)




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
    
    decays = DecayList()
    decayToDM = Decay(instate='Mediator',fstates=['DM','radiation'],br=1.)
    decays.addDecay(decayToDM)
    decays.Xfraction = 0.1
    decays.width = 2.5e-15

    

    #Get the model parameters (or define them here):
    TRH = 1e4
    TF = 1e-3
    
    def dummySigmaV(T,g,mass):
        return g**2*np.exp(-2*mass/T)/(T**2)

    def dummyConvertionRate(T,g,mass,dmass):
        return g**2*np.exp(-dmass/T)*T**2/mass
    
    data = np.genfromtxt('./width_and_medxs.dat',skip_header=5)
    data = np.insert(data,0,[[0.,data[0,1]]],axis=0)
    conv = 0.8579e17
    sigmaVJan = np.vectorize(lambda T: np.exp(interp1d(data[:,0],np.log(data[:,1]*conv),fill_value=0.,bounds_error=False)(500./T)))        

    
    #Define the components to be evolved and their properties:    
    dm = Component(label='DM',Type='thermal',dof=-2,
                   mass=500.,coSigmav=lambda T,other: 1e-12*sigmaVJan(T))
    mediator = Component(label='Mediator',Type='thermal',dof=6,
                   mass=510.,decays=decays,
                   sigmav=sigmaVJan)
    compList = [dm,mediator]
    
    
    #Evolve the equations from TR to TF
    solution = BoltzSolution(compList,TRH,TF,npoints=50000)
    solved = solution.Evolve()
    if not solved:
        return
    
    #Print summary
    if outputFile:
        if os.path.isfile(outputFile):
            os.remove(outputFile)
#         AuxFuncs.printParameters(parser.items('parameters'),outputFile)
        solution.printSummary(outputFile)
        solution.printData(outputFile)
#     else:
    solution.printSummary()
    
    if showPlot:
        #Plot solutions
        import matplotlib.pyplot as plt   
        R = solution.solutionDict['R']
        for comp in compList:
#             n = solution.solutionDict['n_'+comp.label][-1]
            rho = solution.solutionDict['rho_'+comp.label]
            plt.plot(R,rho,label=comp.label)
        plt.plot(R,solution.solutionDict['T'],label='T')
        plt.legend()
        plt.yscale('log')
        plt.xscale('log')
        plt.show()

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
