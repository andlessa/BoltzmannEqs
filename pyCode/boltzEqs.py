#!/usr/bin/env python3

"""

.. module:: boltzEqs
    :synopsis: This module defines the Boltzmann equations 

:synopsis: This module defines the Boltzmann equations
:author: Andre Lessa <lessa.a.p@gmail.com>

"""

from pyCode.AuxFuncs import Hfunc, getTemperature, getPressure
from numpy import exp, log
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
import numpy as np
np.seterr(divide='ignore')
np.seterr(over='ignore')
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


    def __getattr__(self, attr):
        """
        If self does not contain the attribute
        and the attribute has not been defined for the components
        return an array with the values for each component.
        It also applies to methods.

        :param attr: Attribute name

        :return: Array with attribute values
        """

        if not all(hasattr(comp,attr) for comp in self.components):
            raise AttributeError("Components do not have attribute ``%s''" %attr)


        val = getattr(self.components[0],attr)
        if not callable(val):
            return np.array([getattr(br,attr) for br in self.components])

        def call(*args, **kw):
            return np.array([getattr(comp, attr)(*args, **kw) for comp in self.components])
        return call


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

        #Store the number of components:
        nComp = len(self.components)

        #Ni = log(n_i/s_0)
        Ni = y[:nComp]
        #R = rho_i/n_i
        Ri = y[nComp:2*nComp]
        #NS = log(S/S_0)
        NS = y[-1]

        #Get temperature from entropy and scale factor:
        T = getTemperature(x,NS)

        logger.debug('RHS: Computing number and energy densities for %i components' %nComp)
        #Current number densities:
        n = np.exp(Ni)
        #Current energy densities:
        rho = n*Ri
        #Special definition for coherently oscillating fields:
        for i,comp in enumerate(self.components):
            if comp.Type == 'CO':
                rho[i] = comp.mass(T)/Ri[i]

        #Compute equilibrium densities:
        neq = self.nEQ(T)

        #Compute ratio of equilibrium densities
        #(helps with numerical instabilities)
        #rNeq[i,j] = neq[i]/neq[j]
        rNeq = np.array([[compi.rNeq(T,compj) for compj in self.components] for compi in self.components])

        #Dictionary with label:index mapping:
        labelsDict = dict([[comp.label,i] for i,comp in enumerate(self.components)])

        #Compute Hubble factor:
        isActive = self.isActive
        H = Hfunc(T,rho,isActive)
        logger.debug('RHS: Done computing component energy and number densities')

        #Auxiliary weights:
        logger.debug('RHS: Computing weights')
        #Effective equilibrium densities and BRs:
        #NXth[i] = N^{th}_i:
        NXth = self.getNXTh(T,n,rNeq,labelsDict)
        #NXYth[i,j] = N^{th}_{ij}:
        NXYth = np.array([[compi.getNXYTh(T,n,rNeq,labelsDict,compj) for compj in self.components] for compi in self.components])
        #Effective branching ratio (Beff[i,j] = B^{eff}_{ij}:
        Beff = np.array([[compi.getTotalBRTo(T,compj) for compj in self.components] for compi in self.components])
        logger.debug('Done computing weights')

        # Derivative for entropy:
        logger.debug('Computing entropy derivative')     
        dNS = 0.
        for i,comp in enumerate(self.components):
            if not isActive[i]: continue
            dNS += comp.getBRX(T)*comp.width(T)*comp.mass(T)*(n[i]-NXth[i])*exp(3.*x - NS)/(H*T)
        if np.isinf(dNS):
            logger.warning("Infinity found in dNS at T=%1.2g. Will be replaced by a large number" %(T))
            dNS = np.nan_to_num(dNS)

        logger.debug('Done computing entropy derivative')

        #Derivatives for the Ni=log(ni/s0) variables:
        logger.debug('Computing Ni derivatives')
        dN = np.zeros(nComp)
        widths = self.width(T)
        masses = self.mass(T)
        #Expansion term:
        RHS = -3*n
        #Decay term:
        RHS += -widths*masses*n/(H*Ri)
        #Inverse decay term:
        RHS += widths*masses*NXth/(H*Ri) #NXth should be finite if i -> j +..
        #Annihilation term:
        sigmaV = self.getSIGV(T)
        RHS += sigmaV*(neq - n)*(neq + n)/H
        #Contributions from other BSM states:
        for i,compi in enumerate(self.components):
            for j,compj in enumerate(self.components):
                # i + j <-> SM + SM:
                sigVij = compi.getCOSIGV(T,compj)
                if sigVij:
                    RHS[i] += (neq[i]*neq[j]-n[i]*n[j])*sigVij/H #Co-annihilation
                # i+i <-> j+j:
                sigVjj = compi.getSIGVBSM(T,compj)
                if sigVjj:
                    RHS[i] += (rNeq[i,j]*n[j]-n[i])*(rNeq[i,j]*n[j]+n[i])*sigVjj/H #sigVjj*rNeq[i,j]**2 should be finite
                # i+SM <-> j+SM:
                cRate = compi.getConvertionRate(T,compj)
                if cRate:
                    RHS[i] += (rNeq[i,j]*n[j]-n[i])*cRate/H #cRate*rNeq[i,j] should be finite
                # j <-> i +SM:
                RHS[i] += Beff[j,i]*masses[j]*widths[j]*(n[j]-NXYth[j,i])/(H*Ri[j]) #NXYth[j,i] should be finite if j -> i +...

            if not self.isActive[i]:
                if RHS[i] < 0.:
                    continue
                else:
                    logger.warning("Inactive component %s is being injected" %compi.label)
            elif RHS[i]:
                dN[i] = np.float64(RHS[i])/np.float64(n[i])
                if np.isinf(dN[i]):
                    logger.warning("Infinity found at in dN[%s] at T=%1.2g. Will be replaced by a large number" %(comp.label,T))
                    dN[i] = np.nan_to_num(dN[i])


        RHS = np.zeros(nComp)
        dR = np.zeros(nComp)
        #Derivatives for the rho/n variables (only for thermal components):
        for i,comp in enumerate(self.components):
            if not isActive[i] or comp.Type == 'CO':
                continue
            mass = masses[i]
            RHS[i] = -3.*getPressure(mass,rho[i],n[i])  #Cooling term
            for j, compj in enumerate(self.components):
                if not isActive[j]: continue
                if j == i: continue
                massj = masses[j]
                widthj = widths[j]
                #Injection and inverse injection terms:
                RHS[i] += widthj*Beff[j,i]*massj*(1./2. - Ri[i]/Ri[j])*(n[j] - NXYth[j,i])/H #NXth[j,i] should finite if j -> i+..

            if RHS[i]:
                dR[i] = np.float64(RHS[i])/np.float64(n[i])
                if np.isinf(dR[i]):
                    logger.warning("Infinity found in dR[%s] at T=%1.2g. Will be replaced by a large number" %(comp.label,T))
                    dR[i] = np.nan_to_num(dR[i])

        dy = np.hstack((dN,dR,[dNS]))
        logger.debug('T = %1.23g, dNi/dx = %s, dRi/dx = %s, dNS/dx = %s' %(T,str(dN),str(dR),str(dNS)))

        return dy
            
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
                
