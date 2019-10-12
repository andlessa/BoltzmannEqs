#!/usr/bin/env python3

"""

.. module:: boltzSolver
    :synopsis: This module contains the main methods for solving the Boltzmann equations 

:synopsis: This module contains the main methods for solving the Boltzmann equations
:author: Andre Lessa <lessa.a.p@gmail.com>

"""

from pyCode.EqFunctions import gSTAR,gSTARS,Tf,dTfdx
from pyCode.printerTools import printData,printSummary
from numpy import pi,sqrt,log,exp
import numpy as np
from scipy import integrate
import logging
import random, time
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)
random.seed('myseed')


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
                            Ln = np.log(neq*(1+y[icomp])/self.norm[icomp])+100
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
        #For simplicity assume radiation domination for checking thermal
        #equilibrium condition:
        MP = 1.22e19
        H = np.sqrt(8.*np.pi**3*gSTAR(T0)/90.)*T0**2/MP
        sigV = self.getSIGV(T0)
        neq = self.nEQ(T0)
        #Thermal equilibrium condition (if >1 particle is in thermal equilibrium):
        thermalEQ = neq*sigV/H
        for i,comp in enumerate(self.components):
            comp.active = True
            if not hasattr(comp,'n') or not len(comp.n):
                if thermalEQ[i] > 1:
                    comp.n = np.array([neq[i]])
                    comp.coupled = True
                else:
                    comp.n = np.array([1e-20*neq[i]])
                    comp.coupled = False
                comp.rho = comp.n*comp.rEQ(T0)
            else: #Remove all entries from components except the last
                if len(comp.n) > 1:
                    logger.info("Resetting %s's number and energy densities to last value")
                comp.n = comp.n[-1:] 
                comp.rho = comp.rho[-1:]

        #Set thermal equilibrium flags:
        neq = self.nEQ(self.T[-1])
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


            rhoRad = (pi**2/30)*gSTAR(T)*T**4  # thermal bath's energy density    
            MP = 1.22e19
            H = sqrt(8*pi*rhoRad/3)/MP            
            print('sigmav*neq/H=',self.getSIGV(T)*self.nEQ(T)/H)
            
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
        else:
            if r.status < 0:
                logger.error(r.message)
                return False
    
        return r
    
    
    def rhs(self,x,y):
        """
        Defines the derivatives of the y variables at point x = log(R/R0).
        active = [True/False,...] is a list of switches to activate/deactivate components
        If a component is not active it does not evolve and its decay and
        energy density does not contribute to the other components.
        For simplicity we set  R0 = s0 = 1 (with respect to the notes).
        """

        isActive = self.active
        isCoupled = self.coupled
        logger.debug('Calling RHS with arguments:\n   x=%s,\n   y=%s\n and switches %s, %s' %(x,y,isActive,isCoupled))

        #Store the number of components:
        nComp = len(self.components)
        #Auxiliar matrices:
        zeros = np.zeros(nComp)

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
        
        logger.debug('RHS: Computing number and energy densities for %i components' %nComp)
        #Current number densities (replace number densities by equilibrium value for thermally coupled components)
        n = np.array([nv*(1+Ni[i]) if isCoupled[i] else self.norm[i]*np.exp(Ni[i]) 
                      for i,nv in enumerate(neq)])

        #Set number densities to zero if component is not active
        n = np.where(isActive,n,zeros)
        
        #Make sure Ri is never below the component's mass:
        masses = self.mass(T)
        Ri = np.where(Ri > masses,Ri,masses) 
        
        #Define the quantities:
        #   Delta_i = n_i - neq_i, delta_i = 0, for decoupled components
        #   Delta_i = 0, delta_i = n_i/neq_i-1 = y_i, for coupled components
        Di = np.where((not isCoupled)*isActive,n-neq,0.)
        di = np.where(isCoupled*isActive,Ni,0.)
        #Auxilar equilibrium density:
        nT = np.where(isActive,Di+neq,0.)

        #Current energy densities:
        rho = n*Ri

        #Compute ratio of equilibrium densities
        #(helps with numerical instabilities)
        #rNeq[i,j] = neq[i]/neq[j]
        rNeq = np.array([[compi.rNeq(T,compj) if compi.active and compj.active else 0. for compj in self.components] 
                         for compi in self.components])

        #Dictionary with label:index mapping:
        labelsDict = dict([[comp.label,i] for i,comp in enumerate(self.components)])

        #Compute Hubble factor:
        rhoTot = np.sum(rho,where=isActive,initial=0)
        rhoRad = (pi**2/30)*gSTAR(T)*T**4  # thermal bath's energy density    
        rho = rhoRad+rhoTot
        MP = 1.22e19
        H = sqrt(8*pi*rho/3)/MP
        
        logger.debug('RHS: Done computing component energy and number densities')
        logger.debug('n = %s, rho = %s, neq = %s' %(n,rho,neq))

        #Auxiliary weights:
        logger.debug('RHS: Computing process rates and weights')
        
        widths = self.width(T)
        BRX = self.getBRX(T)
        #Process rates:
        sigmaV = self.getSIGV(T) #Annihilation (i+i ->SM+SM)
        sigVij = np.zeros((nComp,nComp)) #Co-annihilation (i+j ->SM+SM)
        sigVjj = np.zeros((nComp,nComp)) #BSM conversion rate (j+j ->i+i)
        cRate = np.zeros((nComp,nComp)) #conversion rate (j+SM ->i+SM)
        #Auxiliary weights:
        NXth = self.getNXTh(T,n,rNeq,labelsDict) #NXth[i] = N^{th}_{i}
        NXYth = np.zeros((nComp,nComp)) #NXYth[i,j] = N^{th}_{ij}
        Beff = np.zeros((nComp,nComp)) #Beff[i,j] = B^{eff}_{ij}       
        for i,compi in enumerate(self.components):
            for j,compj in enumerate(self.components):
                if i == j:
                    continue
                if compi.active and compj.active:
                    sigVij[i,j] = compi.getCOSIGV(T,compj)
                    sigVjj[i,j] = compi.getSIGVBSM(T,compj)
                    cRate[i,j] = compi.getConvertionRate(T,compj)
                if compi.active:
                    NXYth[i,j] = compi.getNXYTh(T,n,rNeq,labelsDict,compj)
                    Beff[i,j] = compi.getTotalBRTo(T,compj)        

        logger.debug('Done computing process rates and weights')

        # Derivative for entropy:
        logger.debug('Computing entropy derivative')     
        dNS = np.sum(isActive*BRX*widths*masses*(n-NXth))*exp(3.*x - NS)/(H*T*self.normS)
        if np.isinf(dNS):
            logger.warning("Infinity found in dNS at T=%1.2g. Will be replaced by a large number" %(T))
            dNS = np.nan_to_num(dNS)

        logger.debug('Done computing entropy derivative')

        #Derivatives for the Ni=log(ni/s0) variables:
        logger.debug('Computing Ni derivatives')
        
        #Terms for decoupled components/zero order terms for coupled components:
        #Expansion term:
        expTerm = -3*(nT*H)
        #Decay term:
        decTerm = -widths*masses*nT/Ri
        #Inverse decay term:
        invDecTerm = widths*masses*NXth/Ri
        #Annihilation term:
        annTerm = Di*(nT+neq)*sigmaV
        coAnnTerm = np.zeros(nComp)
        convTerm = np.zeros(nComp)
        cRateTerm = np.zeros(nComp)
        injectionTerm =np.zeros(nComp)
        for i,compi in enumerate(self.components):
            for j,compj in enumerate(self.components):
                if i == j:
                    continue
                if not compi.active or not compj.active:
                    continue
                # i + j <-> SM + SM:                 
                coAnnTerm[i] += (neq[i]*neq[j]-nT[i]*nT[j])*sigVij[i,j]  
                # i+i <-> j+j (sigVjj*rNeq[i,j]**2 should be finite)
                convTerm[i] += (Di[j]*rNeq[i,j]*(2*neq[i] + Di[j]*rNeq[i,j])
                                - Di[i]*(2*neq[i]+Di[i]))*sigVjj[i,j]
                # i+SM <-> j+SM (cRate*rNeq[i,j] should be finite)
                cRateTerm[i] += (Di[j]*rNeq[i,j] - Di[i])*cRate[i,j]
                # j <-> i +SM (#NXYth[j,i] should be finite if j -> i +...)
                injectionTerm[i] += Beff[j,i]*widths[j]*(masses[j]/Ri[j])*(nT[j] - NXYth[j,i])
        
        #Add all contributions
        allTerms = expTerm+decTerm+invDecTerm+annTerm
        allTerms += coAnnTerm+convTerm+cRateTerm+injectionTerm
        #Make sure non-active components do not evolve:
        allTerms *= isActive
        RHS0 = np.zeros(nComp) 
        np.divide(allTerms,nT*H,where=(allTerms != 0.),out=RHS0)
        #Subtract expansion term for thermally coupled components:
        dLndT = self.dLnEQdT(T) #derivative of log of equilibrium density
        dTdx = dTfdx(x,NS,self.normS)
        RHS0 -= np.where(isCoupled,dTdx*dLndT,zeros)

        logger.debug('DNi (Zero order): expTerm = %1.2g, decTerm = %1.2g, invDecTerm = %1.2g' %(expTerm/(nT*H),decTerm,invDecTerm))
        logger.debug('\t\t annTerm = %1.2g, coAnnTerm = %1.2g, convTerm = %1.2g' %(annTerm,coAnnTerm,convTerm))
        logger.debug('\t\t cRateTerm = %1.2g, injectionTerm = %1.2g, thermalExp = %1.2g' %(cRateTerm,injectionTerm,dTdx*dLndT*isCoupled))
        
        
        #First order terms for coupled components only:
        #Inverse decay term:
        invDecTerm = di*widths*masses*NXth/Ri
        #Annihilation term:
        annTerm = di*2*neq**2*sigmaV
        coAnnTerm = np.zeros(nComp)
        convTerm = np.zeros(nComp)
        cRateTerm = np.zeros(nComp)
        injectionTerm =np.zeros(nComp)
        for i,compi in enumerate(self.components):
            for j,compj in enumerate(self.components):
                if i == j:
                    continue
                if not compi.active or not compj.active:
                    continue
                # i + j <-> SM + SM:                 
                coAnnTerm[i] += -di[i]*neq[i]*neq[j]*sigVij[i,j]
                coAnnTerm[i] += -di[j]*neq[i]*neq[j]*sigVij[i,j]
                # i+i <-> j+j (sigVjj*rNeq[i,j]**2 should be finite)
                convTerm[i] += -di[i]*rNeq[i,j]**2*(nT[j]**2 + neq[j]**2)*sigVjj[i,j]
                convTerm[i] += di[j]*2*neq[i]**3*sigVjj[i,j]/nT[i]
                # i+SM <-> j+SM (cRate*rNeq[i,j] should be finite)
                cRateTerm[i] += -di[i]*nT[j]*rNeq[i,j]*cRate[i,j]
                cRateTerm[i] += di[j]*neq[i]**2*cRate[i,j]/nT[i]
                # j <-> i +SM (#NXYth[j,i] should be finite if j -> i +...)
                injectionTerm[i] += -di[i]*Beff[j,i]*widths[j]*(masses[j]/Ri[j])*(nT[j] - NXYth[j,i])
                injectionTerm[i] += di[j]*Beff[j,i]*widths[j]*(masses[j]/Ri[j])*neq[i]*neq[j]/nT[i]

        logger.debug('DNi (First): invDecTerm = %1.2g' %(invDecTerm))
        logger.debug('\t\t annTerm = %1.2g, coAnnTerm = %1.2g, convTerm = %1.2g' %(annTerm,coAnnTerm,convTerm))
        logger.debug('\t\t cRateTerm = %1.2g, injectionTerm = %1.2g' %(cRateTerm,injectionTerm))

        #Add all contributions
        allTerms = invDecTerm+annTerm
        allTerms += coAnnTerm+convTerm+cRateTerm+injectionTerm
        #Make sure non-active components do not evolve:
        allTerms *= isActive
        RHS1 = np.zeros(nComp) 
        np.divide(allTerms,neq*H,where=(allTerms != 0.),out=RHS1)

        #Add all contributions:
        dN = RHS0+RHS1

        #Derivatives for the rho/n variables (only for thermal components):
        #Compute pressure/n:
        Pn = np.array([comp.Pn(T,Ri[i]) if comp.active else 0. 
                       for i,comp in enumerate(self.components)])
        
        #Terms for decoupled components/zero order terms for coupled components:
        #Expansion/cooling term
        expTerm = -3*(nT*H)*Pn
        injectionTerm =np.zeros(nComp)
        for i,comp in enumerate(self.components):
            for j,compj in enumerate(self.components):            
                if i == j:
                    continue
                if not compi.active or not compj.active:
                    continue
                #Injection and inverse injection terms:
                injectionTerm[i] += Beff[j,i]*widths[j]*masses[j]*(1./2. - Ri[i]/Ri[j])*(nT[j] - NXYth[j,i]) #NXth[j,i] should finite if j -> i+..

        #Add all contributions
        allTerms = expTerm+injectionTerm
        #Make sure non-active components do not evolve:
        allTerms *= isActive
        RHS0 = np.zeros(nComp) 
        np.divide(allTerms,nT*H,where=(allTerms != 0.),out=RHS0)


        #First order terms for coupled components only:
        injectionTerm =np.zeros(nComp)
        for i,comp in enumerate(self.components):
            for j,compj in enumerate(self.components):            
                if i == j:
                    continue
                if not compi.active or not compj.active:
                    continue
                #Injection and inverse injection terms:
                injectionTerm[i] += di[j]*Beff[j,i]*widths[j]*masses[j]*(1./2. - Ri[i]/Ri[j])*neq[j]
                injectionTerm[i] += -di[i]*Beff[j,i]*widths[j]*masses[j]*(1./2. - Ri[i]/Ri[j])*(nT[i]/neq[i])*(nT[j] - NXYth[j,i])

        #Add all contributions
        allTerms = injectionTerm
        #Make sure non-active components do not evolve:
        allTerms *= isActive
        RHS1 = np.zeros(nComp) 
        np.divide(allTerms,nT*H,where=(allTerms != 0.),out=RHS1)


        #Add all contributions:
        dR = RHS0+RHS1

        dy = np.hstack((dN,dR,[dNS]))
        logger.debug('T = %1.23g, dNi/dx = %s, dRi/dx = %s, dNS/dx = %s' %(T,str(dN),str(dR),str(dNS)))

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
        
        #Store the number and energy densities for each component:
        for icomp,comp in enumerate(self.components):
            if not comp.active:
                n = np.array([np.nan]*len(r.t))
                rho = np.array([np.nan]*len(r.t))
            else:
                if not comp.coupled:
                    n = exp(r.y[icomp,:])*self.norm[icomp]
                else:
                    n = (1+r.y[icomp,:])*neq[:,icomp]
                rho = n*r.y[icomp+self.ncomp,:]
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
        
        Delta = y[icomp]
        Delta *= self.components[icomp].active
        Delta *= self.components[icomp].coupled
    
                
        return Delta-0.1
