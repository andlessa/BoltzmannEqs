#!/usr/bin/env python3

#Example to describe the main steps required to define the model and inpute paramaters
#and solve the boltzmann equations

#First tell the system where to find the modules:
import sys,os
from configparser import ConfigParser
import logging
# from matplotlib import pyplot as plt

# logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


def main(parameterFile,outputFile,showPlot=True):
    """
    
    Main code to define the BSM contents and properties and solve the Boltzmann equations
    
    :param parameterFile: Path to the file defining the main model parameters
    :param outputFile: Path to the output file. If None, no results will be written.
    :param showPlot: If True, will show a simple plot for the evolution of the energy densities
    
    """

    import modelDefinitions
    from pyCode.component import Component
    from pyCode.boltzSolver import Evolve

    
    parser = ConfigParser(inline_comment_prefixes=(';',))    
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
    xvals,yvals = Evolve(compList,TRH,TF)
    from matplotlib import pyplot as plt
    plt.plot(xvals,yvals[0,:],'b--',label='DM')
    plt.plot(xvals,yvals[1,:],'r--',label='Mediator')
    plt.show()

   
    return True
    
    
    
    
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
