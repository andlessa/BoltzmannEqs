#!/usr/bin/env python3

"""

.. module:: boltzSolver
    :synopsis: This module contains the main methods for solving the Boltzmann equations 

:synopsis: This module contains the main methods for solving the Boltzmann equations
:author: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys
from pyCode.boltzEqs import BoltzEqs
from pyCode.AuxFuncs import gSTARS, getTemperature, getOmega, getDNeff
from numpy import log,pi,exp
import numpy as np
from scipy import integrate
import logging
import random, time
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)
random.seed('myseed')


class BoltzSolution(object):
    
    def __init__(self,compList,T0,TF,npoints=100,omegaErr=0.01):
        
        self.compList = compList #List of components (species)
        self.T0 = T0 #Initial temperature
        self.TF = TF #Final temperature
        self.omegaErr = omegaErr #Accepted error for relic densities
        self.npoints = npoints
        self.ncomp = len(compList)
        self.solutionDict = {'S' : np.array([]), 'T' : np.array([]), 
                             'x' : np.array([]), 'R' : np.array([])}
        for comp in compList:
            self.solutionDict['n_'+comp.label] = np.array([])
            self.solutionDict['rho_'+comp.label] = np.array([])
        self.boltz_eqs = None
        
    def loadEquations(self):
        """
        Sets the initial conditions and create the BoltzEquations object.
        """
        
        compList = self.compList
        T0 = self.T0
            
        #Sanity checks
        if not self.goodCompList(compList,T0):
            return False
        
        t0 = time.time()
        #Compute initial conditions    
        x0 = 0. # Initial condition for log(R/R0)
        initConditions = np.array([comp.getInitialCond(T0,self.compList) for comp in self.compList])
        y0 = np.concatenate((initConditions[:,0],initConditions[:,1])).tolist()  #Initial conditions for log(n/s0) + conditions for log(rho/n)
        S = (2*pi**2/45)*gSTARS(T0)*T0**3
        y0.append(log(S))  #Initial condition for log(S/S0)
        isActive = [comp.active for comp in compList]
        logger.info("Initial conditions computed in %s s" %(time.time()-t0))
        t0 = time.time()
        
        self.boltz_eqs = BoltzEqs(compList,x0,y0,isActive) #Define equations and set initial conditions
        
    def Evolve(self):
        """
        Evolve the components in component list from the re-heat temperature T0 to TF
        For simplicity we set  R0 = s0 = 1 (with respect to the notes).
        The solution is stored in self.solutionDict.
        Returns True/False if the integration was successful (failed)    
        """
        
        if not self.boltz_eqs:
            self.loadEquations()
        
        t0 = time.time()
        T0 = self.T0
        TF = self.TF
        
        #Solve differential equations:
        x0 = 0.
        #Estimate the final value of x (x=log(R/R0)) assuming conservation of entropy
        xf = log(T0/TF) + (1./3.)*log(gSTARS(T0)/gSTARS(TF)) #Evolve till late times
        tvals = np.linspace(x0,xf,self.npoints)
        y0 = self.boltz_eqs.y0
        logger.debug('Evolving from',x0,'to',xf,'with',len(tvals),'points')
        maxstep = (xf-x0)/300.
        r = integrate.solve_ivp(self.boltz_eqs.rhs,t_span=(x0,xf),y0=y0,
                                t_eval=tvals,method='BDF',
                                events=self.boltz_eqs.events,max_step=maxstep)
        if r.status < 0:
            NS = r.y[-1][-1]
            T = getTemperature(r.t[-1],NS)
            logger.error("Solution failed at temperature %1.3g" %T)
            logger.error("Error message from solver: %s" %r.message)
            return False
        
        self.updateSolution(r)
        x = r.t[-1]
        NS = r.y[-1][-1]
        T  = getTemperature(x,NS)
        while T > TF+maxstep and r.status >= 0:
            x0old = x0            
            xf = x + log(T/TF) + (1./3.)*log(gSTARS(T)/gSTARS(TF)) #Update final x-value for evolutiontvals = np.linspace(x0,xf,self.npoints)            
            if r.t_events:                
                x0,y0,isActive = self.boltz_eqs.handle_events(r)
            else:
                x0 = r.t[-1]
                y0 = r.y[:,-1]
            self.boltz_eqs.updateValues(x0, y0, isActive)
            tvals = np.linspace(x0,xf,int(self.npoints*(xf-x0)/(xf-x0old))+2) #Try to keep the total number of points evenly distributed in x
            logger.debug('Evolving from %1.3g to %1.3g with %i points' %(x0,xf,len(tvals)))
            maxstep = (xf-x0)/300.
            r = integrate.solve_ivp(self.boltz_eqs.rhs,t_span=(x0,xf),y0=y0,
                                        t_eval=tvals,method='BDF',
                                        events=self.boltz_eqs.events,max_step=maxstep)
            self.updateSolution(r)
            x = r.t[-1]
            NS = r.y[-1][-1]
            T  = getTemperature(x,NS)
                          
        
        logger.info("Solution computed in %1.2f s" %(time.time()-t0))
        if r.status < 0:
            logger.error(r.message)
            return False
    
        return True
                    
    def goodCompList(self,compList,T0):
        """Checks if the list of components satisfies the minimum requirements"""
        
        if type(compList) != type(list()):
            logger.error("Input must be a list of Component objects")
            return False
        for comp1 in compList:
            BRs1 = comp1.getBRs(T0)
            for comp2 in compList:            
                if comp1 == comp2: continue
                if not comp2.label in BRs1.getAllFinalStates(): continue
                if comp1.width(T0) < comp2.width(T0):
                    logger.error("Improper decay from %s to %s" %(comp1.label,comp2.label))
                    return False
    
        return True
        
    def updateSolution(self,r):
        """
        Updates the solution in self.solutionDict if the
        integration was successful.
        :param r: Return of scipy.solution_ivp (Bunch object) with information about the integration
        """
        
        if r.status < 0:
            return #Do nothing

        npts = len(r.t)
        
        #Store x-values
        xvalues = r.t
        self.solutionDict['x'] = np.concatenate((self.solutionDict['x'],xvalues))
        
        #Store R values:
        R = exp(r.t)
        self.solutionDict['R'] = np.concatenate((self.solutionDict['R'],R))
        
        #Store T-values
        NSvalues = r.y[-1,:]
        Tvalues = np.array([getTemperature(x, NSvalues[i]) for i,x in enumerate(xvalues)])
        self.solutionDict['T'] = np.concatenate((self.solutionDict['T'],Tvalues))
        
        #Store the number and energy densities for each component:
        for icomp,comp in enumerate(self.compList):
            if not self.boltz_eqs.isActive[icomp]:
                n = np.array([np.nan]*npts)
                rho = np.array([np.nan]*npts)
            else:
                n = exp(r.y[icomp,:])
                rho = n*r.y[icomp+self.ncomp,:]
            self.solutionDict['n_'+comp.label] = np.concatenate((self.solutionDict['n_'+comp.label],n))
            self.solutionDict['rho_'+comp.label] = np.concatenate((self.solutionDict['rho_'+comp.label],rho))
            
        #Store the entropy values:
        S0 = (2*pi**2/45)*gSTARS(self.T0)*self.T0**3
        S = S0*exp(r.y[-1,:])
        self.solutionDict['S'] = np.concatenate((self.solutionDict['S'],S))

    def printSummary(self,outFile=None):
        """
        Prints basic summary of solutions.
        """
        #Solution summary:
        if outFile:
            if hasattr(outFile,'write'):
                f = outFile    
            else:    
                f = open(outFile,'a')
        else:
            f = sys.stdout
            
        T = self.solutionDict['T']
        TF = T[-1]
        f.write('#-------------\n')
        f.write('# Summary\n')
        f.write('# TF=%1.2g\n' %TF)
        for comp in self.compList:
            if comp.Tdecay:
                Tlast = max(comp.Tdecay,TF)
            else:
                Tlast = TF
            #Get point closest to Tlast    
            i = (np.abs(T - Tlast)).argmin()                    
            rhoF = self.solutionDict['rho_'+comp.label][i]
            nF = self.solutionDict['n_'+comp.label][i]
            Tfinal = T[i]
            omega = getOmega(comp,rhoF,nF,Tfinal)        
            if not comp.Tdecay:
                tag = '(@TF)'
            else:
                tag = '(@decay)'
            f.write('# %s: T(osc)= %s | T(decouple)~= %s | T(decay)~= %s | Omega h^2 %s = %1.2f\n' %(comp.label,comp.Tosc,
                                                                                          comp.Tdecouple,comp.Tdecay,tag,omega))
            f.write('# \n')
        
        DNeff = 0.
        for comp in self.compList:
            rho = self.solutionDict['rho_'+comp.label][-1]
            n = self.solutionDict['n_'+comp.label][-1]
            DNeff += getDNeff(comp, rho, n, TF)
                    
        f.write('# Delta Neff (T = %1.2g) = %1.2g \n' %(TF,DNeff))
        f.write('#-------------\n')
        f.close()
    
    def printData(self,outFile=None):
        """
        Prints the evolution of number and energy densities of the species to the outputFile 
        """
        
        if outFile:
            if hasattr(outFile,'write'):
                f = outFile    
            else:    
                f = open(outFile,'a')
            header = sorted(self.solutionDict.keys())
            values = [self.solutionDict[label] for label in header]
            maxLength = max([len(s) for s in header])
            header = ' '.join(str(x).center(maxLength) for x in header)
            if any(len(v) != len(values[0]) for v in values):
                logger.error("Data has distinct lengths and can not be written to file.")
                f.close()
                return False
            data = np.column_stack(values)
            f.write('#--------------\n')
            np.savetxt(f,data,delimiter=' ',header = header,fmt=('{:^%i}'%(maxLength-5)).format('%1.4E'))
            f.write('#--------------\n')    
            f.close()
        
