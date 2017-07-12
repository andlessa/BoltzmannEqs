#!/usr/bin/env python

#Example to describe the main steps required to define the model and inpute paramaters
#and solve the boltzmann equations

#First tell the system where to find the modules:
import sys
from pyCode import AuxFuncs
import modelDefinitions
from pyCode.component import Component
from pyCode.boltzSolver import Evolve
from scipy import integrate, special
from numpy import inf,sqrt,pi,arange

#Fix/change model parameters:
TRH = 1e7
modelDefinitions.yCoupling = (1.5e-12)*1.7
modelDefinitions.mMediator = 500.
modelDefinitions.mDM = 100.
TF = 10.**(-5)


def kint(x,m):
    
    return special.kn(1,x)*x**3/(AuxFuncs.gSTARS(m/x)*sqrt(AuxFuncs.gSTAR(m/x)))

def omegaAnalytic(mDM,mMediator,gMed,width,TRH):    
        
    Ys = (8.488e17)*gMed*width/mMediator**2
    Ttoday = 2.3697*10**(-13)*2.725/2.75  #Temperature today
    Ys *= integrate.quad(kint,mMediator/TRH,inf,mMediator,limit=1000,epsrel=0.01)[0]    
    s = (2.*pi**2*AuxFuncs.gSTARS(Ttoday)*Ttoday**3)/45.
    rhos = mDM*Ys*s    
    omega = rhos/8.0992e-47
    
    om = 1.09e27*gMed*mDM*width/(sqrt(AuxFuncs.gSTAR(mMediator))*AuxFuncs.gSTARS(mMediator)*mMediator**2)
    om *= (AuxFuncs.gSTARS(0.)/3.9091)
    print  'om=',om

    return omega

omegaNum = []
omegaFIMP = []
ypts = []
for x in arange(0.01,5.,0.5):
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
    print 'omegaNum=',omega
    omegaNum.append(omega)
    
    gstar = AuxFuncs.gSTAR(mediator.mass(1.))
    gstars = AuxFuncs.gSTARS(mediator.mass(1.))
    mMediator = mediator.mass(1.)
    mDM = dm.mass(1.)
    
    omega = omegaAnalytic(mDM,mMediator,abs(mediator.dof),mediator.width(1.),TRH)
    
    print 'analytic=',omega
    sys.exit()
    
    omegaFIMP.append(omega)
    ypts.append(modelDefinitions.yCoupling)


print 'FIMP:',omegaFIMP[-1],'Numerical:',omegaNum[-1]
#Plot solutions
import pylab      
pylab.plot(ypts,omegaFIMP,'b',ypts,omegaNum,'r--')
pylab.show()
sys.exit()
