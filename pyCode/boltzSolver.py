#!/usr/bin/env python3

"""

.. module:: boltzSolver
    :synopsis: This module contains the main methods for solving the Boltzmann equations 

:synopsis: This module contains the main methods for solving the Boltzmann equations
:author: Andre Lessa <lessa.a.p@gmail.com>

"""

from pyCode.EqFunctions import gSTAR,gSTARS,Tf
from pyCode.printerTools import printData,printSummary
from numpy import pi,sqrt,log,exp
from scipy import integrate
import numpy as np
import logging
import random, time
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)
random.seed('myseed')
np.set_printoptions(formatter={'float': '{: 1.3g}'.format})

class BoltzSolution(object):
    
    def __init__(self,compList,T0):
        
        self.components = compList #List of components (species)
        self.ncomp = len(compList)
        
        #Define discontinuity events:            
        self.events = []
        for i,comp in enumerate(self.components):
            if not comp.active:
                self.events.append(lambda x,y: -1)
                self.events.append(lambda x,y: -1)
            else:
                def fSuppressed(x,y,icomp=i):
                    comp = self.components[icomp]
                    if comp.active:
                        if comp.coupled:
                            T = Tf(x,y[-1],self.normS)
                            neq = comp.nEQ(T)
                            Ln = np.log(neq/self.norm[icomp])+100
                        else:
                            Ln = y[icomp]+100
                        return Ln
                    else:
                        return 1.
                def fEquilibrium(x,y,icomp=i):
                    return self.checkThermalEQ(x,y,icomp)
                self.events.append(fSuppressed) #Stop evolution particle if its number is too small
                self.events.append(fEquilibrium) #Stop evolution if particle is decoupling
        #Set flag so integration stops when any event occurs:
        for evt in self.events:
            evt.terminal = True
        
        t0 = time.time()
        #Compute initial values for variables:
        self.S = np.array([(2*pi**2/45)*gSTARS(T0)*T0**3])
        self.T = np.array([T0])
        self.x = np.array([0.])
        self.R = np.array([1.])
                
        #Guess initial values for the components (if not yet defined).
        #For simplicity assume radiation domination when checking for thermal
        #equilibrium:
        MP = 1.22e19
        T = self.T[-1]
        H = np.sqrt(8.*np.pi**3*gSTAR(T)/90.)*T**2/MP
        n = neq = self.nEQ(T)
        Ri = self.rEQ(T)
        #Compute all contributions which enforce thermal coupling
        thermalEQ = self.dnidx(T, n, Ri, neq, H, order=1)/neq
        for i,comp in enumerate(self.components):
            comp.active = True
            if not hasattr(comp,'n') or not len(comp.n):
                logger.info("Guessing initial condition for particle %s" %comp)
                if thermalEQ[i] > 1:
                    comp.n = np.array([neq[i]])
                    comp.coupled = True
                else:
                    comp.n = np.array([1e-20*neq[i]])
                    comp.coupled = False
                comp.rho = comp.n*Ri[i]
            else: #Remove all entries from components except the last
                if len(comp.n) > 1:
                    logger.info("Resetting %s's number and energy densities to last value")
                comp.n = comp.n[-1:] 
                comp.rho = comp.rho[-1:]

        #Set thermal equilibrium flags:
        for i,comp in enumerate(self.components):        
            if comp.n[-1] != neq[i] and abs(comp.n[-1]-neq[i])/neq[i] > 1e-4:
                logger.info("Particle %s starting decoupled" %comp)
                comp.coupled = False
            else:
                logger.info("Particle %s starting in thermal equilibrium" %comp)
                comp.coupled = True                
            
        logger.info("Initial conditions computed in %s s" %(time.time()-t0))
        
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

    def setInitialCond(self):
        """
        Use the last entries in the entrory and the components number and energy density values
        to compute the initial conditions.
        """

        Ni0 = np.zeros(self.ncomp) #Initial conditions for Ni = log(ni/ni0) or (ni-neq)/neq for coupled components
        #Initial conditions for Ri = rhoi/ni (set to 1, if component is not active)
        Ri0 = np.where(self.active,self.rho[:,-1]/self.n[:,-1],
                       np.full(self.ncomp,1.))
        NS0 = 0. #Initial condition for NS = log(S/S0)
        self.y0 = np.hstack((Ni0,Ri0,[NS0])).tolist()
        self.norm = self.n[:,-1] #Set normalization (ni0) for each component
        self.normS = self.S[-1] #Set normalization (S0) for entropy
        
    def EvolveTo(self,TF,npoints=5000,dx=None,atol=1e-6,rtol=1e-3):
        """
        Evolve the components in component list from the re-heat temperature T0 to TF
        For simplicity we set  R0 = s0 = 1 (with respect to the notes).
        The solution is stored in self.solutionDict.
        Returns True/False if the integration was successful (failed)    
        """

        #Check if evolution is forward:
        if TF > self.T[-1]:
            logger.error("Final temperature is smaller than initial temperature.")
            return None
        
        #Set initial conditions for equations:
        self.setInitialCond()

        t0 = time.time()
        #Solve differential equations:
        T0 = self.T[0]
        x0 = self.x[-1]
        #Estimate the final value of x (x=log(R/R0)) assuming conservation of entropy
        xf = log(T0/TF) + (1./3.)*log(gSTARS(T0)/gSTARS(TF)) #Evolve till late times
        tvals = np.linspace(x0,xf,npoints)
        y0 = self.y0
        logger.info('Evolving from %1.3g to %1.3g with %i points' %(x0,xf,len(tvals)))
        maxstep = 1e-2
        if dx:
            maxstep = (xf-x0)/dx

        try:
            self.rhs(x0,y0)
        except Exception as e:
            logger.error("Failed to evaluate equations at first point (x0 = %1.3g)." %x0)
            logger.error("Error: %s" %str(e))
            return False
            
        r = integrate.solve_ivp(self.rhs,t_span=(x0,xf),y0=y0,atol=atol,rtol=rtol,
                                t_eval=tvals,method='BDF',dense_output=True,
                                events=self.events,max_step=maxstep)

        if r.status < 0:
            NS = r.y[-1][-1]
            T = Tf(r.t[-1],NS,self.normS)
            logger.error("Solution failed at temperature %1.3g" %T)
            logger.error("Error message from solver: %s" %r.message)
            return False


        self.updateSolution(r)
        logger.info("Solution computed in %1.2f s" %(time.time()-t0))        
        
        continueEvolution = False        
        for i,evt in enumerate(r.t_events):
            comp = self.components[int(i/2)]
            if evt.size > 0:
                continueEvolution = True
                if np.mod(i,2):
                    logger.info("Integration restarted because %s left thermal equilibrium at x=%s (T = %1.3g GeV)" %(comp.label,str(evt),self.T[-1]))
                    comp.Tdecouple = self.T[-1]
                    comp.coupled = False
                else:
                    logger.info("Integration restarted because the number density for %s became too small at x=%s (T = %1.3g GeV)" %(comp.label,str(evt),self.T[-1]))
                    comp.Tdecay = self.T[-1]
                    comp.active = False
        if continueEvolution and any(comp.active for comp in self.components):
            self.EvolveTo(TF, npoints-len(r.t), dx, atol, rtol)
        elif all(not comp.active for comp in self.components):
            logger.info("Evolution stopping since all components are inactive.")
        
        if r.status < 0:
            logger.error(r.message)
            return False
    
        return r
    
    def getVariables(self,x,y):
        
        isActive = self.active
        isCoupled = self.coupled

        #Store the number of components:
        nComp = len(self.components)

        #Ni = log(n_i/s_0)
        Ni = np.array(y[:nComp])
        #R = rho_i/n_i
        Ri = np.array(y[nComp:2*nComp])
        #NS = log(S/S_0)
        NS = y[-1]

        #Get temperature from entropy and scale factor:
        T = Tf(x,NS,self.normS)

        #Compute equilibrium densities:
        neq = self.nEQ(T)
        Req = self.rEQ(T)
        
        #Current number densities (replace number densities by equilibrium value for thermally coupled components)
        n = np.where(isCoupled,neq,self.norm*np.exp(Ni))
        #Set number densities to zero if component is not active
        n = np.where(isActive,n,0.)
        
        #Set energy density ratio to equilibrium value for thermally coupled components:
        Ri = np.where(isCoupled,Req,Ri)
        #Make sure Ri is never below the component's mass:
        masses = self.mass(T)
        Ri = np.maximum(Ri,masses)
        
        variables = np.array([T,n,Ri,NS,neq,Req])

        return variables
    
    def getRates(self,T):
        
        #Store the number of components:
        nComp = len(self.components)
                
        #Process rates:
        sigVii = self.sigVii(T) #Annihilation (i+i ->SM+SM)
        sigVij = np.zeros((nComp,nComp)) #Co-annihilation (i+j ->SM+SM)
        sigVjj = np.zeros((nComp,nComp)) #BSM conversion rate (j+j ->i+i)
        cRate = np.zeros((nComp,nComp)) #conversion rate (j+SM ->i+SM)
        
        for i,compi in enumerate(self.components):
            for j,compj in enumerate(self.components):
                if i == j:
                    continue
                if compi.active and compj.active:
                    sigVij[i,j] = compi.sigVij(T,compj)
                    sigVjj[i,j] = compi.sigVjj(T,compj)
                    cRate[i,j] = compi.cRate(T,compj)
       
        return sigVii,sigVij,sigVjj,cRate
    
    def getAuxNumbers(self,T,n,rNeq,labelsDict):
        
        #Store the number of components:
        nComp = len(self.components)
        
        #Auxiliary weights:
        NXth = self.getNXTh(T,n,rNeq,labelsDict) #NXth[i] = N^{th}_{i}
        NXYth = np.zeros((nComp,nComp)) #NXYth[i,j] = N^{th}_{ij}
        Beff = np.zeros((nComp,nComp)) #Beff[i,j] = B^{eff}_{ij}       
        for i,compi in enumerate(self.components):
            for j,compj in enumerate(self.components):
                if i == j:
                    continue
                if compi.active:
                    NXYth[i,j] = compi.getNXYTh(T,n,rNeq,labelsDict,compj)
                    Beff[i,j] = compi.getTotalBRTo(T,compj)        
        
        return NXth,NXYth,Beff
    
    def dNSdx(self,T,x,n,NS,H):
        """
        Compute derivatives for the log of total entropy.
        """

        isActive = self.active
        rNeq = np.array([[compi.rNeq(T,compj) if compi.active and compj.active else 0. 
                          for compj in self.components] for compi in self.components])
        labelsDict = dict([[comp.label,i] for i,comp in enumerate(self.components)])
        
        BRX = self.getBRX(T)
        masses = self.mass(T)
        widths = self.width(T)
        NXth,_,_ = self.getAuxNumbers(T,n,rNeq,labelsDict)        
        
        # Derivative for entropy:
        dNS = np.sum(isActive*BRX*widths*masses*(n-NXth))*exp(3.*x - NS)/(H*T*self.normS)
        if np.isinf(dNS):
            logger.warning("Infinity found in dNS at T=%1.2g. Will be replaced by a large number" %(T))
            dNS = np.nan_to_num(dNS)
            
        return dNS

    def dnidx(self,T,n,Ri,neq,H,order=0):
        """
        Compute derivatives for the number density.
        """
        
        isActive = self.active
        #Store the number of components:
        nComp = len(self.components)
        rNeq = np.array([[compi.rNeq(T,compj) if compi.active and compj.active else 0. 
                          for compj in self.components] for compi in self.components])
        labelsDict = dict([[comp.label,i] for i,comp in enumerate(self.components)])  

        masses = self.mass(T)
        widths = self.width(T)
        sigVii,sigVij,sigVjj,cRate = self.getRates(T)
        NXth,NXYth,Beff =  self.getAuxNumbers(T,n,rNeq,labelsDict)
        isActiveM = np.outer(isActive,isActive)
        diagM = np.identity(nComp)
        offM = np.ones((nComp,nComp))-diagM
        neqMatrix = np.outer(neq,neq)
        nMatrix = np.outer(n,n)

        if order == 0: #Compute all terms contributing to dni/dx:
            #Expansion term:
            expTerm = -3*n
            #Decay term:
            decTerm = -widths*masses*n/(H*Ri)
            #Inverse decay term:
            invDecTerm = widths*masses*NXth/(H*Ri)
            #Annihilation term (i +i <-> SM + SM):
            annTerm = np.where(sigVii*isActive,sigVii*(neq - n)*(neq + n)/H,0.)
            #Co-annihilation term (i + j <-> SM + SM):
            coAnnTerm = np.where(sigVij*isActiveM*offM,(neqMatrix-nMatrix)*sigVij,0.)
            coAnnTerm = np.sum(coAnnTerm,axis=1)/H
            #Conversion annilihation (i+i <-> j+j) (sigVjj*rNeq[i,j]**2 should be finite)
            convTerm = np.where(sigVjj*isActiveM*offM,(rNeq**2*n**2-n[:,np.newaxis]**2)*sigVjj,0.)
            convTerm = np.sum(convTerm,axis=1)/H
            #Conversion rate (i+SM <-> j+SM) (cRate*rNeq[i,j] should be finite)
            cRateTerm = np.where(cRate*isActiveM*offM,(rNeq*n-n[:,np.newaxis])*cRate,0.) 
            cRateTerm = np.sum(cRateTerm,axis=1)/H
            #Injection term (j <-> i +SM) (#NXYth[j,i] should be finite if j -> i +...)
            injectionTerm = np.where(Beff.T*widths*isActiveM*offM,Beff.T*widths*(masses/Ri)*(n-NXYth.T),0.)
            injectionTerm = np.sum(injectionTerm,axis=1)/H
            
            #Add all contributions
            allTerms = expTerm+decTerm+invDecTerm+annTerm
            allTerms += coAnnTerm+convTerm+cRateTerm+injectionTerm
            
        elif order == 1: #Compute all terms in dnidx which are first order in delta_i = (n/neq)-1
            expTerm = decTerm = 0.
            annTerm = 2*neq**2*sigVii/H
            invDecTerm = widths*masses*NXth/(H*Ri)
            coAnnTerm =  neq*sigVij.dot(neq)/H
            convTerm = (sigVjj*rNeq**2).dot(n**2 + neq**2)/H
            cRateTerm = (cRate*rNeq).dot(n)/H
            injectionTerm = np.where(Beff.T*widths*isActiveM*offM,Beff.T*widths*(masses/Ri)*(n-NXYth.T),0.)
            injectionTerm = np.sum(injectionTerm,axis=1)/H
            
            #Add all contributions which enforce thermal equilibrium
            allTerms = invDecTerm+annTerm
            allTerms += coAnnTerm+convTerm+cRateTerm+injectionTerm
        
        else:
            logger.error("Order must either be zero or one (not %s)" %order)
            return None

        logger.debug('DNi: expTerm = %s, decTerm = %s, invDecTerm = %s' %(expTerm,decTerm,invDecTerm))
        logger.debug('\t\t annTerm = %s, coAnnTerm = %s, convTerm = %s' %(annTerm,coAnnTerm,convTerm))
        logger.debug('\t\t cRateTerm = %s, injectionTerm = %s' %(cRateTerm,injectionTerm))
           
        return allTerms
 
    def dRidx(self,T,n,Ri,H):
        """
        Compute derivatives for the ratio of energy and number densities.
        """
        
        isActive = self.active
        #Store the number of components:
        nComp = len(self.components)
        rNeq = np.array([[compi.rNeq(T,compj) if compi.active and compj.active else 0. 
                          for compj in self.components] for compi in self.components])
        labelsDict = dict([[comp.label,i] for i,comp in enumerate(self.components)])  

        masses = self.mass(T)
        widths = self.width(T)
        _,NXYth,Beff =  self.getAuxNumbers(T,n,rNeq,labelsDict)
        isActiveM = np.outer(isActive,isActive)
        diagM = np.identity(nComp)
        offM = np.ones((nComp,nComp))-diagM
        
        #Compute pressure/n:
        Pn = np.array([comp.Pn(T,Ri[i]) if comp.active else 0. 
                       for i,comp in enumerate(self.components)])
        
        #Expansion/cooling term
        expTerm = -3*Pn*n
        #Approximate energy injection (matrix):
        eInjectionM = (1./2.)*(1+np.outer(masses**2,1/masses**2))
        rMatrix = np.outer(Ri,1/Ri)
        injectionTerm = np.where(Beff.T*widths*isActiveM*offM,
                                 Beff.T*widths*masses*(eInjectionM-rMatrix)*(n-NXYth.T),
                                 0.)
        injectionTerm = np.sum(injectionTerm,axis=1)/H
        #Add all contributions
        allTerms = expTerm+injectionTerm
        #Make sure non-active and coupled components do not evolve:
        dR = np.where(allTerms*isActive,allTerms/n,0.)
                
        logger.debug('DRi: expTerm = %s, injectionTerm = %s' %(expTerm,injectionTerm))
        
        return dR

    def rhs(self,x,y):
        """
        Get the derivatives of the y variables at point x = log(R/R0).
        active = [True/False,...] is a list of switches to activate/deactivate components
        If a component is not active it does not evolve and its decay and
        energy density does not contribute to the other components.
        For simplicity we set  R0 = s0 = 1 (with respect to the notes).
        """

        isActive = self.active
        isCoupled = self.coupled
        logger.debug('Calling RHS with arguments:\n   x=%s,\n   y=%s\n and isActive = %s, isCoupled = %s' %(x,y,isActive,isCoupled))

        #Store the number of components:
        nComp = len(self.components)
        
        logger.debug('RHS: Computing number and energy densities for %i components' %nComp)
        
        T,n,Ri,NS,neq,Req = self.getVariables(x, y)

        #Compute Hubble factor:
        rho = n*Ri
        rhoTot = np.sum(rho,where=isActive,initial=0)
        rhoRad = (pi**2/30)*gSTAR(T)*T**4  # thermal bath's energy density    
        rhoTot += rhoRad
        MP = 1.22e19
        H = sqrt(8*pi*rhoTot/3)/MP
        
        logger.debug('n = %s, rho = %s, neq = %s, Req = %s' %(n,rho,neq,Req))
        
        # Derivative for entropy:
        logger.debug('Computing entropy derivative')     
        dNS = self.dNSdx(T, x, n, NS, H)

        #Derivatives for the Ni=log(ni/s0) variables:
        logger.debug('Computing Ni derivatives')
        allTerms = self.dnidx(T, n, Ri, neq, H)
        #Make sure non-active and coupled components do not evolve:
        dN = np.where(allTerms*isActive*np.invert(isCoupled),allTerms/n,0.)

        #Derivatives for the rho/n variables:
        logger.debug('Computing Ri derivatives')
        dR = self.dRidx(T, n, Ri, H)
        #Make sure non-active and coupled components do not evolve:
        dR *= isActive*np.invert(isCoupled)
        
        dy = np.hstack((dN,dR,[dNS]))
        logger.debug('T = %1.3g, H = %1.3g, dNi/dx = %s, dRi/dx = %s, dNS/dx = %s' %(T,H,dN,dR,dNS))

        return dy

    def updateSolution(self,r):
        """
        Updates the solution in self.solutionDict if the
        integration was successful.
        :param r: Return of scipy.solution_ivp (Bunch object) with information about the integration
        """
        
        if r.status < 0:
            return #Do nothing

        #Store x-values
        self.x = np.hstack((self.x,r.t))       
        #Store R values:
        self.R = np.hstack((self.R,np.exp(r.t)))
        #Store the entropy values:
        S = self.normS*exp(r.y[-1,:])
        self.S = np.hstack((self.S,S))  
        #Store T-values
        NSvalues = r.y[-1,:]
        Tvalues = np.array([Tf(x,NSvalues[i],self.normS) for i,x in enumerate(r.t)])
        self.T = np.hstack((self.T,Tvalues))
        neq = np.array([self.nEQ(T) for T in Tvalues])
        Req = np.array([self.rEQ(T) for T in Tvalues])
        
        #Store the number and energy densities for each component:
        for icomp,comp in enumerate(self.components):
            if not comp.active:
                n = np.array([np.nan]*len(r.t))
                rho = np.array([np.nan]*len(r.t))
            else:
                if not comp.coupled:
                    n = exp(r.y[icomp,:])*self.norm[icomp]
                    rho = n*r.y[icomp+self.ncomp,:]
                else:
                    n = neq[:,icomp]
                    rho = n*Req[:,icomp]
                
            comp.n = np.hstack((comp.n,n))
            comp.rho = np.hstack((comp.rho,rho))

    def printSummary(self,outFile=None):
        """
        Prints basic summary of solutions.
        """
        printSummary(self, outFile)
    
    def printData(self,outFile=None):
        """
        Prints the evolution of number and energy densities of the species to the outputFile 
        """
        
        printData(self, outFile)

    def checkThermalEQ(self,x,y,icomp):
        """
        Check if component icomp is leaving thermal equilibrium.
        """
        
        isActive = self.active
        isCoupled = self.coupled
        
        if not isActive[icomp]*isCoupled[icomp]:
            return 1.0
        
        #Store the number of components:
        T,n,Ri,_,neq,_ = self.getVariables(x, y)
        
        #Compute Hubble factor:
        rho = n*Ri
        rhoTot = np.sum(rho,where=isActive,initial=0)
        rhoRad = (pi**2/30)*gSTAR(T)*T**4  # thermal bath's energy density    
        rhoTot += rhoRad
        MP = 1.22e19
        H = sqrt(8*pi*rhoTot/3)/MP
        

        #Compute all contributions which enforce thermal coupling
        allTermsEq = self.dnidx(T,n,Ri,neq,H,order=1)[icomp]
        if not allTermsEq:
            return 0.
        
        #Compute all contributions which might lead to thermal decoupling
        allTermsNonEq = self.dnidx(T,n,Ri,neq,H,order=0)[icomp]

        #Consider as softly coupled when non-equilibrium contributions 
        #are 1% of the thermal equilibrium ones
        r = 1. - 100*abs(allTermsNonEq/allTermsEq)
        
        return r
    
