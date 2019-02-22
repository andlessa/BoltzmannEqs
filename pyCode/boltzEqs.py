#!/usr/bin/env python

"""

.. module:: boltzEqs
    :synopsis: This module defines the Boltzmann equations 

:synopsis: This module defines the Boltzmann equations
:author: Andre Lessa <lessa.a.p@gmail.com>

"""

from pyCode.AuxFuncs import Hfunc, getTemperature, getPressure
from math import isnan, pi
import logging
from scipy.integrate import ode
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
import numpy as np
import decimal
from scipy.special import zetac
np.seterr(over='ignore')
Zeta3 = zetac(3.) + 1.
decimal.getcontext().prec = 28


class BoltzEqs(object):
    """Class to hold the Boltzmann equations"""
    
#Sets the initial conditions and check for errors in input
    def __init__(self,compList,x0=None,y0=None,sw=None):
        self.components = compList
        self.solver = ode(self.rhs)
        self.solver.set_integrator('lsoda',method='bdf')
        self.sw = sw
        if not x0 is None and not y0 is None:
            self.solver.set_initial_value(y0, x0)


    #The right-hand-side function (rhs)
    def rhs(self,x,y):
        """
        Defines the derivatives of the y variables at point x = log(R/R0).
        sw = [True/False,...] is a list of switches to activate/deactivate components
        If a component is not active it does not evolve and its decay and
        energy density does not contribute to the other components.
        For simplicity we set  R0 = s0 = 1 (with respect to the notes).
        """

#         self.check_discontinuities(x,y)
        sw = self.sw
        logger.debug('Calling RHS with arguments:\n   x=%s,\n   y=%s\n and switches %s' %(x,y,sw))
        n = []
        neq = []
        rho = []
        R = []
        nratio = {'radiation' : 1.}   #n/neq ratio dictionary
        NS = y[-1]
        T = getTemperature(x,NS)
        for i,comp in enumerate(self.components):
            logger.debug('RHS: Computing component %s' %comp)
            ni = np.exp(y[i])
            Ri = y[i + len(self.components)]
            if comp.Type == 'CO':
                rhoi = comp.mass(T)*ni
            else: rhoi = Ri*ni
            n.append(ni)
            neq.append(comp.nEQ(T))            
            rho.append(rhoi)
            R.append(Ri)
            if neq[-1] > 0.: nratio[comp.label] = ni/neq[-1]
            else: nratio[comp.label] = 0.
            logger.debug('RHS: Done computing component %s.\n   rho = %s and n = %s' %(comp,rhoi,ni))
            if y[i] < -100.:
                sw[i] = False
            else:
                sw[i] = True
        H = Hfunc(T,rho,sw)
       
             
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
                if a == i or not sw[a]: continue
                N2th[a][i] = compA.getNTh(T,nratio,comp)
                Beff[a][i] = compA.getTotalBRTo(T,comp)
        logger.debug('Done computing weights')
# Derivative for entropy:
        logger.debug('Computing entropy derivative')     
        dNS = 0.        
        for i,comp in enumerate(self.components):
            if not sw[i]: continue
            dNS += comp.getBRX(T)*comp.width(T)*comp.mass(T)*(n[i]-N1th[i])*np.exp(3.*x - NS)/(H*T)
        logger.debug('Done computing entropy derivative')
             
#Derivatives for the Ni=log(ni/s0) variables:
        logger.debug('Computing Ni derivatives')
        dN = [0.]*len(self.components)        
        for i,comp in enumerate(self.components):
            if not sw[i]:
                continue
            width = comp.width(T)
            mass = comp.mass(T)
            RHS = -3.*n[i]     
            RHS += -width*mass*(n[i] - N1th[i])/(H*R[i])    #Decay term
            RHS += comp.getSource(T)/H  #Source term
            annTerm = 0.
            if comp.Type == 'weakthermal':                
                nrel = Zeta3*T**3/pi**2
                annTerm = comp.getSIGV(T)*nrel/H                
            else:
                annTerm = comp.getSIGV(T)*n[i]/H
            #Define approximate decoupling temperature (just used for printout)            
            if annTerm < 1e-2 and not comp.Tdecouple:
                comp.Tdecouple = T
            elif annTerm > 1e-2 and comp.Tdecouple:
                comp.Tdecouple = None  #Reset decoupling temperature if component becomes coupled
            RHS += annTerm*(neq[i] - n[i]) #Annihilation term
            for a, compA in enumerate(self.components):
                if not sw[a]:
                    continue
                if a == i:
                    continue                                
                massA = compA.mass(T)
                widthA = compA.width(T)
                RHS += widthA*Beff[a][i]*massA*(n[a] - N2th[a][i])/(H*R[a])  #Injection term                       
            dN[i] = RHS/n[i]    #Log equations
        

            

        dR = [0.]*len(self.components)
#Derivatives for the rho/n variables (only for thermal components):        
        for i,comp in enumerate(self.components):
            if not sw[i] or comp.Type == 'CO': continue                 
            mass = comp.mass(T)
            RHS = -3.*getPressure(mass,rho[i],n[i])/n[i]  #Cooling term
            for a, compA in enumerate(self.components):
                if not sw[a]: continue                
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
        
        logger.debug('dN= %s' %dN)
        logger.debug('dR= %s' %dR)
        logger.debug('dNS= %s' %dNS)
        logger.debug('DONE')
        return np.array(dN + dR + [dNS])
    

    def check_discontinuities(self,x,y):
        """
        Activate/de-activate components when a discontinuous transition happens
        and sets the new initial conditions.          
        Possible discontinuities are: a CO component starts to oscillate, a particle has decayed.
        """
        
        #Event happened because of discontinuity
        NS = y[-1]
        T = getTemperature(x,NS)         
        for icomp,yv in enumerate(y):            
            if self.components[icomp].width(T) == 0. or not self.sw[icomp]:
                continue
            if 2.*yv + 200. < 0:
                self.sw[icomp] = False
