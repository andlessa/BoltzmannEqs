#!/usr/bin/env python

#Example to describe the main steps required to define the model and inpute paramaters
#and solve the boltzmann equations

#First tell the system where to find the modules:
import sys
from pyCode import AuxFuncs
import modelDefinitions
from pyCode.component import Component
from pyCode.boltzSolver import Evolve
from numpy import arange
from modelDefinitions import mMediator

#Fix/change model parameters:
TRH = 1e7
modelDefinitions.yCoupling = (1.5e-10)
modelDefinitions.mMediator = 1000.
modelDefinitions.mDM = 1.
TF = 10.**(-8)



def main():
    omegaNum = []
    omegaFIMP = []
    ypts = []
    Tdecay = []
    Tdecouple = []
    npts = 0
    for x in arange(-4.,1.5,0.1):
        modelDefinitions.yCoupling = (1.5e-10)*10**x
        #Define the components to be evolved:   
        dm = Component(label='DM',Type='thermal',dof=1,
                       mass=modelDefinitions.mDM,sigmav=modelDefinitions.DMSigmaV)
        mediator = Component(label='Mediator',Type='thermal',dof=-2,
                       mass=modelDefinitions.mMediator,decays=modelDefinitions.MediatorDecays,
                       sigmav=modelDefinitions.MediatorSigmaV)
        compList = [dm,mediator]    
        #Evolve the equations from TR to TF
        Evolve(compList,TRH,TF)    
        rhoF = dm.evolveVars['rho'][-1]
        nF = dm.evolveVars['n'][-1]
        Tfinal = dm.evolveVars['T'][-1]        
        omega = AuxFuncs.getOmega(dm,rhoF,nF,Tfinal)
        omegaNum.append(omega)

        outputFile = open('./note/yscan/scan_%s.dat' %npts, 'w')
        npts += 1
        parsDict = {'mMediator' : modelDefinitions.mMediator,'mDM' : modelDefinitions.mDM,
                    'TRH' : TRH, 'TF' : TF, 'yCoupling' : modelDefinitions.yCoupling}
        AuxFuncs.printParameters(parsDict.items(),outputFile)
        AuxFuncs.printSummary(compList,TF,outputFile)
#         AuxFuncs.printData(compList,outputFile)
                
        
        mMediator = mediator.mass(1.)
        mDM = dm.mass(1.)
        Tdecay.append(mediator.Tdecay)    
        Tdecouple.append(mediator.Tdecouple)        
        omega = AuxFuncs.omegaAnalytic(mDM,mMediator,abs(mediator.dof),mediator.width(1.),TRH)
        outputFile.write('Omega h^2 (analytic) = %s\n' %omega)
        outputFile.close()
        omegaFIMP.append(omega)
        ypts.append(modelDefinitions.yCoupling)
    
    print 'FIMP:',omegaFIMP[-1],'Numerical:',omegaNum[-1]
#     #Plot solutions
#     import pylab
#     omdiff = [omegaNum[i]/om for i,om in enumerate(omegaFIMP)]
#     # pylab.plot(ypts,omegaFIMP,'b',ypts,omegaNum,'r--')
#     Tratio = [T/Tdecay[i] for i,T in enumerate(Tdecouple)]
#     pylab.loglog(ypts,omdiff,'r-o',ypts,Tratio,'b-o')
#     pylab.show()
#     
    return

if __name__ == "__main__":
    main()

