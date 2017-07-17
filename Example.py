#!/usr/bin/env python

#Example to describe the main steps required to define the model and inpute paramaters
#and solve the boltzmann equations

#First tell the system where to find the modules:
import sys,os
from ConfigParser import SafeConfigParser
import logging as logger


def main(parameterFile,outputFile,showPlot=True):
    """
    
    Main code to define the BSM contents and properties and solve the Boltzmann equations
    
    :param parameterFile: Path to the file defining the main model parameters
    :param outputFile: Path to the output file. If None, no results will be written.
    :param showPlot: If True, will show a simple plot for the evolution of the energy densities
    
    """

    from pyCode import AuxFuncs
    import modelDefinitions
    from pyCode.component import Component
    from pyCode.boltzSolver import Evolve

    
    parser = SafeConfigParser()    
    if not parser.read(parameterFile):
        logger.error("No such file or directory: '%s'" % parameterFile)
        sys.exit()

    #Get the model parameters (or define them here):
    TRH = parser.getfloat("parameters","TRH")
    TF = parser.getfloat("parameters","TF")
    modelDefinitions.yCoupling = parser.getfloat("parameters","yCoupling")
    modelDefinitions.mMediator = parser.getfloat("parameters","mMediator")
    modelDefinitions.mDM = parser.getfloat("parameters","mDM")
    
    #Define the components to be evolved and their properties:   
    dm = Component(label='DM',Type='thermal',dof=1,
                   mass=modelDefinitions.mDM,sigmav=modelDefinitions.DMSigmaV)
    mediator = Component(label='Mediator',Type='thermal',dof=-2,
                   mass=modelDefinitions.mMediator,decays=modelDefinitions.MediatorDecays,
                   sigmav=modelDefinitions.MediatorSigmaV)
    compList = [dm,mediator]
    
    #Evolve the equations from TR to TF
    Evolve(compList,TRH,TF)
    
    #Print summary
    # TF = compList[0].evolveVars['T'][-1]    
    if outputFile:
        if os.path.isfile(outputFile):
            os.remove(outputFile)
        AuxFuncs.printParameters(parser.items('parameters'),outputFile)
        AuxFuncs.printSummary(compList,TF,outputFile)
        AuxFuncs.printData(compList,outputFile)
    else:
        AuxFuncs.printSummary(compList,TF,sys.stdout)
    
    if showPlot:
        #Plot solutions
        import pylab      
        for comp in compList:
            pylab.plot(comp.evolveVars['R'],comp.evolveVars['rho'],label=comp.label)
        pylab.plot(compList[0].evolveVars['R'],compList[0].evolveVars['T'],label='T')
        pylab.legend()
        pylab.yscale('log')
        pylab.xscale('log')
        pylab.show()

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
