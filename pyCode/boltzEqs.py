#!/usr/bin/env python3

"""

.. module:: boltzEqs
    :synopsis: This module defines the Boltzmann equations 

:synopsis: This module defines the Boltzmann equations
:author: Andre Lessa <lessa.a.p@gmail.com>

"""

from pyCode.AuxFuncs import Hfunc, getTemperature, getPressure
from numpy import exp, log, isnan
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
import numpy as np
from scipy.special import zetac
Zeta3 = zetac(3.) + 1.
# np.seterr(over='ignore')


class BoltzEqs(object):
    """Class to hold the Boltzmann equations"""
    
    #Sets the initial conditions
    def __init__(self,compList,x0=None,y0=None,isActive=None):
        self.components = compList    
        #Set initial conditions  
        self.updateValues(x0,y0,isActive)
        #Define discontinuity events:            
        self.events = [self.check_decayOf(i) for i in range(len(compList))]
        self.events += [self.check_oscillationOf(i) for i in range(len(compList))]

    def updateValues(self,x,y,isActive):
        """Replace the variables in self.components and set the new initial conditions."""        
        #Store new values:
        if len(y) != 2*len(self.components)+1:
            logger.error("Wrong number of variables. Must be equal to number of 2*components+1 (entropy)")
            return False
        #Set initial conditions for next evolution
        self.t0 = x                              
        self.y0 = y     
        self.isActive = isActive
    
    def rhs(self,x,y):
        """
        Defines the derivatives of the y variables at point x = log(R/R0).
        isActive = [True/False,...] is a list of switches to activate/deactivate components
        If a component is not active it does not evolve and its decay and
        energy density does not contribute to the other components.
        For simplicity we set  R0 = s0 = 1 (with respect to the notes).
        """

        isActive = self.isActive
        logger.debug('Calling RHS with arguments:\n   x=%s,\n   y=%s\n and switches %s' %(x,y,isActive))
        n = []
        neq = []
        rho = []
        R = []
        nratio = {'radiation' : 1.}   #n/neq ratio dictionary
        NS = y[-1]
        T = getTemperature(x,NS)
        for i,comp in enumerate(self.components):
            logger.debug('RHS: Computing component %s' %comp)
            ni = exp(y[i])
            Ri = y[i + len(self.components)]
            if comp.Type == 'CO':
                rhoi = comp.mass(T)*ni
            else:
                rhoi = Ri*ni
            n.append(ni)
            neq.append(comp.nEQ(T))            
            rho.append(rhoi)
            R.append(Ri)
            if neq[-1] > 0.:
                nratio[comp.label] = ni/neq[-1]
            else:
                nratio[comp.label] = 0.
            logger.debug('RHS: Done computing component %s.\n   rho = %s and n = %s' %(comp,rhoi,ni))
        H = Hfunc(T,rho,isActive)
       
             
#Auxiliary weights:
        logger.debug('Computing weights')     
        N1th = [0.]*len(self.components)
        N2th = [0.]*len(self.components)
        Beff = [0.]*len(self.components)
        for i,comp in enumerate(self.components):
            N2th[i] = [0.]*len(self.components)
            Beff[i] = [0.]*len(self.components)                        
        for i,comp in enumerate(self.components):
            N1th[i] = comp.getNTh(T,nratio)
            for a,compA in enumerate(self.components):                                
                if a == i or not isActive[a]: continue
                N2th[a][i] = compA.getNTh(T,nratio,comp)
                Beff[a][i] = compA.getTotalBRTo(T,comp)
        logger.debug('Done computing weights')
# Derivative for entropy:
        logger.debug('Computing entropy derivative')     
        dNS = 0.        
        for i,comp in enumerate(self.components):
            if not isActive[i]: continue
            dNS += comp.getBRX(T)*comp.width(T)*comp.mass(T)*(n[i]-N1th[i])*exp(3.*x - NS)/(H*T)
        logger.debug('Done computing entropy derivative')
             
#Derivatives for the Ni=log(ni/s0) variables:
        logger.debug('Computing Ni derivatives')
        dN = [0.]*len(self.components)        
        for i,comp in enumerate(self.components):
            if not isActive[i]: continue
            width = comp.width(T)
            mass = comp.mass(T)
            RHS = -3.*n[i]
            RHS += -width*mass*(n[i] - N1th[i])/(H*R[i])    #Decay term
            RHS += comp.getSource(T)/H  #Source term
            annTerm = 0. #Annihilation term
            coannTerm = 0. #Co-annihilation term
            bsmScatter = 0. #2<->2 scattering between BSM components
            convertion = 0. #i<->j convertion
            annTerm = comp.getSIGV(T)*(neq[i]**2 - n[i]**2)/H
            for j, compj in enumerate(self.components):
                if j == i or not isActive[j]:
                    continue
                rj = n[j]/neq[j]
                coannTerm += comp.getCOSIGV(T,compj)*(neq[i]*neq[j]-n[i]*n[j])/H
                bsmScatter += comp.getSIGVBSM(T,compj)*(neq[i]**2*rj**2-n[i]**2)/H
                convertion += comp.getConvertionRate(T,compj)*(neq[i]*rj-n[i])/H
            annTotal = annTerm+coannTerm+bsmScatter+convertion 
            #Define approximate decoupling temperature (just used for printout)            
            if annTotal < 1e-10 and not comp.Tdecouple:
                comp.Tdecouple = T
            elif annTotal > 1e-10 and comp.Tdecouple:
                comp.Tdecouple = None  #Reset decoupling temperature if component becomes coupled
            RHS += annTotal #Annihilation term
            for a, compA in enumerate(self.components):
                if not isActive[a]: continue
                if a == i: continue                                
                massA = compA.mass(T)
                widthA = compA.width(T)
                RHS += widthA*Beff[a][i]*massA*(n[a] - N2th[a][i])/(H*R[a])  #Injection term

            dN[i] = RHS/n[i]    #Log equations

            

        dR = [0.]*len(self.components)
#Derivatives for the rho/n variables (only for thermal components):        
        for i,comp in enumerate(self.components):
            if not isActive[i] or comp.Type == 'CO': continue                 
            mass = comp.mass(T)
            RHS = -3.*getPressure(mass,rho[i],n[i])/n[i]  #Cooling term
            for a, compA in enumerate(self.components):
                if not isActive[a]: continue                
                if a == i: continue                
                massA = compA.mass(T)
                widthA = compA.width(T)
                RHS += widthA*Beff[a][i]*massA*(1./2. - R[i]/R[a])*(n[a]/n[i] - N2th[a][i]/n[i])/H  #Injection term
            
            dR[i] = RHS

        bigerror = False
        for val in y:
            if isnan(val) or abs(val) == float('inf'): bigerror = 1     
        for val in np.array(dN + dR + [dNS]):
            if not bigerror and (isnan(val) or abs(val) == float('inf')): bigerror = 2
        if bigerror:
            if bigerror == 1:  logger.warning("Right-hand called with NaN values.")
            if bigerror == 2:  logger.warning("Right-hand side evaluated to NaN.")
        
        print(dN,dR,dNS)
        return np.array(dN + dR + [dNS])
            
    def check_decayOf(self,icomp):
        """
        Create event for checking if component icomp has decayed
        """
        def event(x,y):
            NS = y[-1]
            T = getTemperature(x,NS)
            comp = self.components[icomp]                     
            if comp.width(T) == 0. or not self.isActive[icomp]:
                return 1            
            return y[icomp] + 100.
        
        event.terminal = True #Whether to stop the evolution when event happens
        return event
    
    def check_oscillationOf(self,icomp):
        """
        Create event for checking if component icomp started to oscillate
        """
        
        def event(x,y):
            comp = self.components[icomp]
            if comp.Type != 'CO':
                return 1            
            NS = y[-1]
            T = getTemperature(x,NS)           
            
            ncomp = len(self.components)
            n = exp(np.array(y)[:ncomp])
            rho = n*np.array(y[ncomp:2*ncomp])
            return Hfunc(T,rho,self.isActive)*3. - comp.mass(T)
        
        event.terminal = True #Whether to stop the evolution when event happens
        return event        

    def handle_events(self,solution):
        """Activate/de-activate components when a discontinuous transition happens
        and sets the new initial conditions.
        Possible discontinuities are: a CO component starts to oscillate, a particle has decayed                 
        """
        
#Event happened because of discontinuity
        NS = solution.y[:,-1][-1]
        T = getTemperature(solution.t[-1],NS)
        events = solution.t_events
        ncomp = len(self.components)
        newY0 = solution.y[:,-1]
        newX0 = solution.t[-1]
        isActive = self.isActive[:]
        for icomp,comp in enumerate(self.components):
            decayT, oscillT = None,None
            if events[icomp].size:
                decayT = events[icomp][0]
                T = getTemperature(decayT,NS)
            if events[icomp+ncomp].size:
                oscillT = events[icomp+ncomp][0]
            if not decayT and not oscillT:
                continue
            if oscillT:
                if comp.Type != 'CO':
                    logger.error("Particle started of oscillate, but it is of type %s" %comp.Type)
                    return False
                isActive[icomp] = True
                comp.active = True 
                comp.Tosc = T
                comp.xosc = solution.t[-1]
                initN = comp.getOscAmplitute(T)/comp.mass(T)
                newY0[icomp] = log(initN)  #Initial condition for log(n/s0)
                newY0[icomp+ncomp] = comp.mass(T) #Initial condition for log(rho/n)
            if decayT:    #Discontinuities due to component decay                
                if self.components[icomp].width(T) == 0.:
                    logger.error("Stable component %s has decayed. Should not happen." %self.components[icomp].label)
                    return False
                comp.Tdecay = T
                isActive[icomp] = False   #Deactivate particle
                
        return newX0,newY0,isActive
                
