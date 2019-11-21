#!/usr/bin/env python3

"""

.. module:: component
    :synopsis: This module defines the component class, which describe common properties for particle/field components 
        
:author: Andre Lessa <lessa.a.p@gmail.com>

"""

from pyCode.AuxDecays import DecayList
from pyCode.EqFunctions import Tf,gSTAR,gSTARS,Pnf,rEQf,nEQf,rNeqf,dLnEQfdT,drEQdTf
import numpy as np
from scipy import integrate
from types import FunctionType
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

Types = ['thermal','CO', 'weakthermal']

class Component(object):
    """Main class to hold component properties.
    
    :param label: label for the component (can be anything)
    :param ID: a proper name for defining the component (must belong to IDs)
    :param Type: must be weakthermal, thermal or CO (weakly coupled thermal, thermal or coherent oscillating field)
    :param active: if component is active or not (True/False)
    :param dof: +/- Number of degrees of freedom (must be negative for fermions and positive for bosons)
    :param Tdecay: stores the temperature at which the particle starts to decay (approximate)
    :param Tosc: stores the temperature at which the particle starts to oscillate (if Type = CO)
    :param Tdecouple: stores the temperature at which the particle decouples (if Type = thermal/weakthermal)
    :param evolveVars: Stores the values of the variables used for evolution
    
    """
    
    def __init__(self,label,Type,dof,mass,decays=DecayList(),
                 sigVii = lambda x: 0., sigVij = lambda x,y: 0.,
                 cRate = lambda x,y: 0., sigVjj = lambda x,y: 0.):
        self.label = label
        self.Type = Type
        self.active = True
        self.dof = dof
        self.Tdecay = None
        self.Tosc = None
        self.Tdecouple = None
        self.sigVii = sigVii
        self.sigVij = sigVij
        self.cRate = cRate
        self.sigVjj = sigVjj

        if not Type or type(Type) != type(str()) or not Type in Types:
            logger.error("Please define proper particle Type (not "+str(Type)+"). \n Possible Types are: "+str(Types))
            return False

        #Check if given mass already is a function: 
        if isinstance(mass,FunctionType):
            self.mass = mass #Use function given
        elif isinstance(mass,int) or isinstance(mass,float):
            self.mass = lambda T: float(mass)  #Use value given for all T
        else:
            logger.error("Mass must be a number or a function of T")
            return False

        #Check if given decays already is a function: 
        if isinstance(decays,FunctionType):
            self.decays = decays #Use function given
        elif isinstance(decays,DecayList):
            self.decays = lambda T: decays #Use value given for all T
        else:
            logger.error("Decays must be a DecayList object or a function of T")
            return False

    def __str__(self):
        
        return self.label


    def getBRX(self,T):
        """
        Returns the fraction of energy injected in the thermal bath by the component decay at temperature T
        """ 
       
        return self.decays(T).Xfraction
    
    def width(self,T):
        """
        Returns the component width at temperature T.
        """
          
        return self.decays(T).width
    
    def getTotalBRTo(self,T,comp):
        """
        Computes the total branching ratio to compoenent comp (sum BR(self -> comp + X..)*comp_multiplicity),
        including the multiplicity factor
        """

        if self is comp:
            return 0.

        brTot = 0.
        for decay in self.decays(T):
            if not comp.label in decay.fstateIDs: continue
            brTot += decay.fstateIDs.count(comp.label)*decay.br
        return brTot
    
    def getNXTh(self,T,n,rNeq,labelsDict):
        """        
        Computes the effective thermal number density of first type at temperature T:
          
        :param T: temperature (allows for T-dependent BRs)
        :param n: list of number densities
        :param rNeq: list with ratios of number densities
        :param labelsDict: Dictionary with the component label -> index in n,rNeq mapping

        :return: Effective thermal number density at temperature T of first type (N_X^{th})
        """

        #Compute BRs:
        Nth = 0.
        BRs = self.decays(T)
        #Index for self:
        i = labelsDict[self.label]
        neq = self.nEQ(T)
        #If the thermal equilibrium density of self is zero,
        #there is no inverse decay:
        if not neq:
            return 0.
        #Catch overflow if neq << 1 (in this case inverse decay should be zero):
        try:
            norm = 1/neq
        except RuntimeWarning:
            return 0.
        for decay in BRs:
            nprod = 1.
            norm = 1./neq
            if not decay.br:
                continue  #Ignore decays with zero BRs
            for label in decay.fstateIDs:
                if label in labelsDict:
                    j = labelsDict[label]
                    norm *= neq
                    nprod *= rNeq[i,j]*n[j]
            if not norm:
                return 0.
            Nth += nprod*decay.br/norm

        return Nth

    def getNXYTh(self,T,n,rNeq,labelsDict,comp):
        """
        Computes the effective thermal number density of second type at temperature T
        for self -> comp +X decays.

        :param T: temperature (allows for T-dependent BRs)
        :param n: list of number densities
        :param rNeq: list with ratios of number densities
        :param labelsDict: Dictionary with the component label -> index in n,rNeq mapping
        :param comp: Component object.

        :return: Effective thermal number density at temperature T of second type for self -> comp (N_{XY}^{th})
        """

        Beff = self.getTotalBRTo(T,comp)
        #If self does not decay to comp, return zero
        if not Beff:
            return 0.

        #Compute BRs:
        Nth = 0.
        BRs = self.decays(T)
        #Index for self:
        i = labelsDict[self.label]
        neq = self.nEQ(T)
        #If the thermal equilibrium density of comp is zero,
        #there is no inverse injection:
        if not neq:
            return 0.
        try:
            norm = 1/neq
        except RuntimeWarning:
            return 0.

        for decay in BRs:
            nprod = 1.
            norm = 1/neq
            if not decay.br:
                continue  #Ignore decays with zero BRs
            if not comp.label in decay.fstateIDs:
                continue #Ignore decays which do not include self
            for label in decay.fstateIDs:
                if label in labelsDict:
                    j = labelsDict[label]
                    nprod *= rNeq[i,j]*n[j]
                    norm *= neq

            if not norm:
                return 0.

            Nth += nprod*decay.br*decay.fstateIDs.count(comp.label)/norm

        Nth *= 1./Beff

        return Nth
    
    
    def Pn(self,T,R):
        """
        Computes the ratio of pressure and number density
        for a component given the ratio of its energy and number densities.
        """
        
        mass = self.mass(T)
        
        return Pnf(R,mass)
            
    def nEQ(self,T):
        """
        Returns the equilibrium number density at temperature T.
        """
        
        mass = self.mass(T)
        dof = self.dof
        
        return nEQf(T,mass,dof)  
    
    
    def rNeq(self,T,other):
        """
        Returns the ratio of equilibrium number densities at temperature T nEQ_self/nEQ_other.
        If both species are Boltzmann suppressed (T << mass) and nearly degenerate,
        it provides a better approximation than taking the ratio of equilibrium densities,
        since the Boltzmann suppression exponentials partially cancel.
        """

        if self is other:
            return 1.
        
        massA = self.mass(T)
        massB = other.mass(T)
        dofA = self.dof
        dofB = other.dof
        
        return rNeqf(T,massA,dofA,massB,dofB)

    def rEQ(self,T):
        """
        Returns the ratio of equilibrium energy and number densities at temperature T,
        assuming chemical potential = 0.
        """

        mass = self.mass(T)
        dof = self.dof
    
        return rEQf(T,mass,dof)
    
    def dLnEQdT(self,T):
        """
        Returns the derivative of the equilibrium number density at temperature T
        with respect to T divided by the equilibrium number density.
        """
        
        mass = self.mass(T)
        
        return dLnEQfdT(T,mass)  
    
    def drEQdT(self,T):
        """
        Returns the derivative of the ratio of equilibrium energy and number densities 
        with respect to the temperature T. Assumes zero chemical potential.
        """
        
        mass = self.mass(T)
        dof = self.dof
        
        return drEQdTf(T,mass,dof)  

    def guessInitialCond(self,T,components=[]):
        """
        Get initial conditions for component at temperature T, assuming the energy
        density is dominated by radiation.
        """
        
        MP = 1.22e19
        
        H = np.sqrt(8.*np.pi**3*gSTAR(T)/90.)*T**2/MP
        
        coannTerm = 0. #Co-annihilation term
        bsmScatter = 0. #2<->2 scattering between BSM components
        convertion = 0. #i<->j convertion
        annTerm = self.getSIGV(T)*self.nEQ(T)/H #Annihilation term
        for comp in components:
            if comp is self:
                continue
            if comp.Tdecouple:
                continue
            coannTerm += self.getCOSIGV(T,comp)*self.nEQ(T)/H
            bsmScatter += self.getSIGVBSM(T,comp)*self.nEQ(T)/H
            convertion += self.getConvertionRate(T,comp)/H

        totalProdRate = annTerm+convertion+bsmScatter+coannTerm
        if totalProdRate < 2.:  #Particle starts decoupled
            logger.info("Particle %s starting decoupled" %self)
            return 1e-10*self.nEQ(T)
        else: #Particle starts coupled
            logger.info("Particle %s starting in thermal equilibrium" %self)
            return self.nEQ(T)
           
    def getOmega(self,rho,n,T):
        """
        Compute relic density today, given the component, number density and energy density
        at temperature T.
        """
        
        if self.Tdecay and self.Tdecay > T:
            return 0.
        if not n or not rho:
            return 0.
        
        Ttoday = 2.3697e-13*2.725/2.75  #Temperature today
        rhoh2 = 8.0992e-47   # value of rho critic divided by h^2
        xToday = (1./3.)*np.log(gSTARS(T)/gSTARS(Ttoday)) + np.log(T/Ttoday)   #dx = log(R/R_today), where R is the scale factor
        nToday = n*np.exp(-3.*xToday)
        s0 = (2*np.pi**2/45)*T**3
        ns = 0.  #log(entropy) (constant) 
        
        R0 = rho/n
        Rmin = R0*np.exp(-xToday)    #Minimum value for rho/n(Ttoday) (happens if component is relativistic today)
        Pmin = nToday*Pnf(Ttoday,Rmin)
        
        if abs(Pmin - Rmin*nToday/3.)/(Rmin*nToday/3.) < 0.01:
            RToday = Rmin  #Relativistic component today
        else:
            def Rfunc(R,x):
                TF = Tf(x,ns,s0)                        
                return -3*self.Pn(TF,R)
            RToday = integrate.odeint(Rfunc, R0, [0.,xToday], atol = self.mass(Ttoday)/10.)[1][0]  #Solve decoupled ODE for R=rho/n
       
        return RToday*nToday/rhoh2
    
    def getDNeff(self,rho,n,T):
        """
        Computes the contribution from the component to the number of effective neutrinos at temperature T.
        Gives zero if T > 1 MeV (where neutrinos are still coupled).
        """
    
    #Ignore component if it has decayed before TF        
        if self.Tdecay and self.Tdecay > T:
            return 0.
    #Get the number and energy densities of comp at T:
        mass = self.mass(T)
        if T > 10.**(-3):
            return 0.    
        if mass == 0. or (n and rho and rho/(n*mass) > 2.):
            rhoRel = rho    
        else:
            rhoRel = 0.
        DNeff = rhoRel/(((np.pi**2)/15)*(7./8.)*((4./11.)**(4./3.))*T**4)
        
        return DNeff    
        
