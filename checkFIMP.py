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

#Fix/change model parameters:
TRH = 1e7
modelDefinitions.yCoupling = (1.5e-12)*1.7
modelDefinitions.mMediator = 500.
modelDefinitions.mDM = 100.
TF = 10.**(-8)



def main():
    omegaNum = []
    omegaFIMP = []
    ypts = []
    Tdecay = []
    Tdecouple = []
    for x in arange(0.05,5.,0.5):
        modelDefinitions.yCoupling = (1.5e-12)*x
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
        
        mMediator = mediator.mass(1.)
        mDM = dm.mass(1.)
        Tdecay.append(mediator.Tdecay)    
        Tdecouple.append(mediator.Tdecouple)
        
        omega = AuxFuncs.omegaAnalytic(mDM,mMediator,abs(mediator.dof),mediator.width(1.),TRH)
        omegaFIMP.append(omega)
        ypts.append(modelDefinitions.yCoupling)
    
    
    print 'FIMP:',omegaFIMP[-1],'Numerical:',omegaNum[-1]
    #Plot solutions
    import pylab
    omdiff = [omegaNum[i]/om for i,om in enumerate(omegaFIMP)]
    # pylab.plot(ypts,omegaFIMP,'b',ypts,omegaNum,'r--')
    Tratio = [T/Tdecay[i] for i,T in enumerate(Tdecouple)]
    pylab.loglog(ypts,omdiff,'r-o',ypts,Tratio,'b-o')
    pylab.show()
    
    return

if __name__ == "__main__":
    main()

