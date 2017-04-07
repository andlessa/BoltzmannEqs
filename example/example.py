#!/usr/bin/env python

#Example to describe the main steps required to define the model and inpute paramaters
#and solve the boltzmann equations

#First tell the system where to find the modules:
codedir = "../pyCode/"
import sys,os
sys.path.append(os.path.abspath(codedir))
from math import sqrt
import modelParameters,AuxFuncs

#Data files
sigmav_datafile = "sigmaV.dat"
rge_datafile = "rge.dat"
slhafile = "spectrum.slha"
isajet_inputfile = 'omegasusy.in'
#Load the relevant data from datafiles and run isajer/isared (if needed)
modelParameters.getIsajetData(codedir,sigmav_datafile,slhafile,rge_datafile,isajet_inputfile)

#PQMSSM parameters:
fa = 10.**10
TR = 10.**6
modelParameters.codeDir = codedir
modelParameters.fa = fa
modelParameters.vPQ = fa/sqrt(2.)
modelParameters.thetaI = 1.
modelParameters.sI = fa
modelParameters.cH = 1.
modelParameters.xi = 1.
modelParameters.mass_axino = 1000.
modelParameters.mass_saxion = 1000.
modelParameters.mass_gravitino = 1010.
modelParameters.TR = TR
modelParameters.modelType = 'DFSZ'  #DFSZ or KSVZ
modelParameters.nW = 1  #Number of Heavy quark representations (only for KSVZ)
modelParameters.caYY = 8./3.  #Number of Heavy quark representations (only for KSVZ)



from component import Component
from boltzSolver import Evolve
#Define the components to be evolved:   
axion = Component(label='Axion',ID='axion',Type='thermal',dof=1)
saxion = Component(label='Saxion',ID='saxion',Type='thermal',dof=1)
axionCO = Component(label='Axion (CO)',ID='axionCO',Type='CO',dof=1)
saxionCO = Component(label='Saxion (CO)',ID='saxionCO',Type='CO',dof=1)
axino = Component(label='Axino',ID='axino',Type='thermal',dof=-2)
neutralino = Component(label='Z1',ID='neutralino',Type='thermal',dof=-2)
gravitino = Component(label='Gravitino',ID='gravitino',Type='thermal',dof=-4)
compList = [neutralino,axion,axionCO,saxion,saxionCO,axino,gravitino]


#Evolve the equations from TR to TF
TF = 10.**(-5)
Evolve(compList,TR,TF)

#Print summary
# TF = compList[0].evolveVars['T'][-1]
AuxFuncs.printSummary(compList,TF)

#Plot solutions
import pylab      
for comp in compList:
    pylab.plot(comp.evolveVars['R'],comp.evolveVars['rho'],label=comp.label)
pylab.plot(compList[0].evolveVars['R'],compList[0].evolveVars['T'],label='T')
pylab.legend()
pylab.yscale('log')
pylab.xscale('log')
pylab.show()
sys.exit()


#Save solutions to file
out = open('rho.dat','w')
header = "T:R:rad"
for comp in compList: header += ":"+comp.label
out.write(header+"\n")
Tpts = compList[0].evolveVars['T']  #Common to all components
Rpts =  compList[0].evolveVars['R']
for ipt,T in Tpts:
    toprint = str(T)+" "+str(Rpts[ipt])+" "+str((modelParameters.Pi**2/30.)*AuxFuncs.gSTAR(T)*T**4)  #Temperature,Scale Factor, rho_radiation
    for ic,comp in enumerate(compList):
        n = comp.evolveVars['n'][ipt]
        rho = comp.evolveVars['rho'][ipt]
        toprint += " "+str(rho)  #Print only energy densities
    out.write(toprint+"\n")       
out.close()
