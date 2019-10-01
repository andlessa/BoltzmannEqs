#!/usr/bin/env python3

"""

.. module:: boltzSolver
    :synopsis: This module contains the main methods for solving the Boltzmann equations 

:synopsis: This module contains the main methods for solving the Boltzmann equations
:author: Andre Lessa <lessa.a.p@gmail.com>

"""

from pyCode.EqFunctions import gSTAR, gSTARf,gSTARSf,Tf
from pyCode.EqFunctions import T as Tfunc
from pyCode.printerTools import printData,printSummary
import sympy as sp
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
                    if self.components[icomp].active:
                        return y[icomp]+100.
                    else:
                        return 1.
                def fEquilibrium(x,y,icomp=i):
                    return 1.
#                     return self.checkThermalEQ(x,y,icomp)
                self.events.append(fSuppressed) #Stop evolution particle if its number is too small
                self.events.append(fEquilibrium) #Stop evolution if particle is decoupling
        #Set flag so integration stops when any event occurs:
        for evt in self.events:
            evt.terminal = True
        
        t0 = time.time()
        #Compute initial values for variables:
        self.S = np.array([(2*np.pi**2/45)*gSTARSf(T0)*T0**3])
        self.T = np.array([T0])
        self.x = np.array([0.])
        self.R = np.array([1.])
                
        #Guess initial values for the components (if not yet defined).
        #For simplicity assume radiatio            n domination for checking thermal
        #equilibrium condition:
        MP = 1.22*10**19
        H = np.sqrt(8.*np.pi**3*gSTARf(T0)/90.)*T0**2/MP
        sigV = self.sigmavF(T0)
        neq = self.nEQf(T0)
        #Thermal equilibrium condition (if >1 particle is in thermal equilibrium):
        thermalEQ = neq*sigV/H
        for i,comp in enumerate(self.components):
            comp.active = True
            if not hasattr(comp,'n') or not len(comp.n):
                if thermalEQ[i] > 1:
                    comp.n = np.array([neq[i]])
                else:
                    comp.n = np.array([1e-20*neq[i]])
                    comp.Tdecouple = T0
                comp.rho = comp.n*comp.rEQf(T0)
            else: #Remove all entries from components except the last
                if len(comp.n) > 1:
                    logger.info("Resetting %s's number and energy densities to last value")
                comp.n = comp.n[-1:]             
                comp.rho = comp.rho[-1:]

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

        Ni0 = np.zeros(self.ncomp) #Initial conditions for Ni = log(ni/ni0)
        Ri0 = self.rho[:,-1]/self.n[:,-1] #Initial conditions for Ri = rhoi/ni
        NS0 = 0. #Initial condition for NS = log(S/S0)
        self.y0 = np.hstack((Ni0,Ri0,[NS0])).tolist()
        self.norm = self.n[:,-1] #Set normalization (ni0) for each component
        self.normS = self.S[-1] #Set normalization (S0) for entropy
        
    def EvolveTo(self,TF,npoints=5000,dx=None,doJacobian=True,
                 atol=1e-6,rtol=1e-3):
        """
        Evolve the components in component list from the re-heat temperature T0 to TF
        For simplicity we set  R0 = s0 = 1 (with respect to the notes).
        The solution is stored in self.solutionDict.
        Returns True/False if the integration was successful (failed)    
        """
        
        #Set initial conditions for equations:
        self.setInitialCond()
        
        #Define Boltzmann equations and the jacobian (if required)
        self.rhs, self.jac = self.getRHS(doJacobian)

        t0 = time.time()
        #Solve differential equations:
        T0 = self.T[0]
        x0 = self.x[-1]
        #Estimate the final value of x (x=log(R/R0)) assuming conservation of entropy
        xf = np.log(T0/TF) + (1./3.)*np.log(gSTARSf(T0)/gSTARSf(TF)) #Evolve till late times
        tvals = np.linspace(x0,xf,npoints)
        y0 = self.y0
        logger.debug('Evolving from %1.3g to %1.3g with %i points' %(x0,xf,len(tvals)))
        maxstep = np.inf
        if dx:
            maxstep = (xf-x0)/dx
        r = integrate.solve_ivp(self.rhs,t_span=(x0,xf),y0=y0,atol=atol,rtol=rtol,
                                t_eval=tvals,method='BDF',dense_output=True,
                                events=self.events,max_step=maxstep,jac=self.jac)
        
        if r.status < 0:
            NS = r.y[-1][-1]
            Tfail = Tf(r.t[-1],NS,self.normS)
            logger.error("Solution failed at temperature %1.3g" %Tfail)
            logger.error("Error message from solver: %s" %r.message)
            return False


        self.updateSolution(r)
        
        continueEvolution = False        
        for i,evt in enumerate(r.t_events):
            comp = self.components[int(i/2)]
            if evt.size > 0:
                continueEvolution = True
                if np.mod(i,2):
                    logger.info("Integration restarted because %s left thermal equilibrium at x=%s (T = %1.3g GeV)" %(comp.label,str(evt),self.T[-1]))
                    comp.Tdecouple = self.T[-1]
                else:
                    logger.info("Integration restarted because the number density for %s became too small at x=%s (T = %1.3g GeV)" %(comp.label,str(evt),self.T[-1]))
                    comp.Tdecay = self.T[-1]
                    comp.active = False
        
        if continueEvolution and any(comp.active for comp in self.components):
            self.EvolveTo(TF, npoints-len(r.t), dx)
        else:
            logger.info("Solution computed in %1.2f s" %(time.time()-t0))                          
            if r.status < 0:
                logger.error(r.message)
                return False
    
        return r
    
    def getRHS(self,doJacobian=True):
        
        """
        Obtain numerical functions for evaluating the right-hand side
        of the Boltzmann equations and its Jacobian (if required).
        First the algebraic equations are defined using Ni,Ri,NS and
        x as sympy variables. Derivatives of discontinuos or non-analytic
        functions are not explicitly computed and only evaluated numerically.
        
        :param doJacobian: Boolean specifying if the Jacobian of the differential
                           equations should be computed. If False, will return
                           None for the jacobian function.
                           
        :return: The RHS and Jacobian functions to be evaluated numerically.
        
        """
        
        #Store the number of components:
        nComp = len(self.components)

        #Define variables:        
        N = np.array(sp.symbols('N:%d'%nComp))
        R = np.array(sp.symbols('R:%d'%nComp))
        NS = sp.symbols('N_S')
        x = sp.symbols('x')
        T = Tfunc(x,NS,self.normS)
        
        #Planck constant:
        MP = sp.symbols('M_P')
        
        #Current number densities:
        n = self.norm*np.array([sp.exp(Ni) for Ni in N])
        #Current energy densities:
        rho = n*R
        
        #Compute equilibrium densities:
        neq = self.nEQ(T)
        
        #Compute ratio of equilibrium densities
        #(helps with numerical instabilities)
        #rNeq[i,j] = neq[i]/neq[j]
        rNeq = np.array([[compi.rNeq(T,compj) if compi.active and compj.active else 0. for compj in self.components] 
                         for compi in self.components])
        
        #Dictionary with label:index mapping:
        labelsDict = dict([[comp.label,i] for i,comp in enumerate(self.components)])
        isActive = self.active
        
        #Compute Hubble factor:
        rhoTot = np.sum(rho,where=isActive,initial=0)
        rhoRad = (sp.pi**2/30)*gSTAR(T)*T**4  # thermal bath's energy density    
        rho = rhoRad+rhoTot
        H = sp.sqrt(8*sp.pi*rho/3)/MP
                
        #Auxiliary weights:
        #Effective equilibrium densities and BRs:
        #NXth[i] = N^{th}_i:
        NXth = self.getNXTh(T,n,rNeq,labelsDict)
        #NXYth[i,j] = N^{th}_{ij}:
        NXYth = np.array([[compi.getNXYTh(T,n,rNeq,labelsDict,compj) if not compj is compi and compj.active else 0. for compj in self.components] 
                          for compi in self.components])
        #Effective branching ratio (Beff[i,j] = B^{eff}_{ij}:
        Beff = np.array([[compi.getTotalBRTo(T,compj) if not compj is compi and compj.active else 0. for compj in self.components] 
                         for compi in self.components])        
        widths = self.width(T)
        masses = self.mass(T)
        BRX = self.getBRX(T)
        sigmaV = self.getSIGV(T)
        
        # Derivative for entropy:
        dNS = np.sum(isActive*BRX*widths*masses*(n-NXth))*sp.exp(3.*x - NS)/(H*T*self.normS)
        
        #Derivatives for the Ni=log(ni/s0) variables:
        #Expansion term:
        RHS = -3*n
        #Annihilation term:            
        RHS += sigmaV*(neq - n)*(neq + n)/H
        # i + j <-> SM + SM:
        for i,compi in enumerate(self.components):
            for j,compj in enumerate(self.components):
                if i == j or not compj.active:
                    continue
                sigVij = compi.getCOSIGV(T,compj)
                if not sigVij:
                    continue
                RHS += (neq[i]*neq[j]-n[i]*n[j])*sigVij/H
        # i+i <-> j+j (sigVjj*rNeq[i,j]**2 should be finite)
        for i,compi in enumerate(self.components):
            for j,compj in enumerate(self.components):
                if i == j or not compj.active:
                    continue
                sigVjj = compi.getSIGVBSM(T,compj)
                if not sigVjj:
                    continue
                RHS += (rNeq[i,j]**2*n[j]**2-n[i]**2)*sigVjj/H
        # i+SM <-> j+SM (cRate*rNeq[i,j] should be finite)
        for i,compi in enumerate(self.components):
            for j,compj in enumerate(self.components):
                if i == j or not compj.active:
                    continue
                cRate = compi.getConvertionRate(T,compj)
                if not cRate:
                    continue
                RHS += (rNeq[i,j]*n[j]-n[i])*cRate/H
        # j <-> i +SM (#NXYth[j,i] should be finite if j -> i +...)
        for i,compi in enumerate(self.components):
            for j,compj in enumerate(self.components):
                if i == j or not compj.active:
                    continue
                cRate = compi.getConvertionRate(T,compj)
                if not cRate:
                    continue
                RHS += Beff[j,i]*widths[j]*masses[j]*(n[j]-NXYth[j,i])/(H*R[j])
        #Decay and inverse decay terms:
        RHS -= widths*masses*(n-NXth)/(H*R) #NXth should be finite if i -> j +.. 


        for i,compi in enumerate(self.components):
            if not isActive[i] and RHS[i] > 0.:
                logger.warning("Inactive component %s is being injected" %compi.label)
        
        dN = sp.sympify(np.zeros(nComp)).as_mutable()
        for i,rhs in enumerate(RHS):
            if isActive[i]:
                dN[i] = rhs/n[i]


        #Derivatives for the Ri=rhoi/ni variables:
        RHS = sp.sympify(np.zeros(nComp)).as_mutable()
        #Derivatives for the rho/n variables (only for thermal components):
        for i,comp in enumerate(self.components):
            if not isActive[i]:
                continue
            RHS[i] = -3.*n[i]*comp.Pn(T,R[i])  #Cooling term
            for j, compj in enumerate(self.components):
                if i == j or not compj.active:
                    continue
                #Injection and inverse injection terms:
                RHS[i] += Beff[j,i]*widths[j]*masses[j]*(1./2. - R[i]/R[j])*(n[j] - NXYth[j,i])/H #NXth[j,i] should finite if j -> i+..

        dR = sp.sympify(np.zeros(nComp)).as_mutable()
        for i,rhs in enumerate(RHS):
            if isActive[i]:
                dR[i] = rhs/n[i]

        dy = np.hstack((dN,dR,[dNS])) #Derivatives
        yv = np.hstack((N,R,[NS])) #y-variables
        
        #Convert the algebraic equation in a numerical equation:
        rhsf = sp.lambdify([x,yv],dy, 
                          modules=[{'M_P' : 1.22e19},labelsDict,'numpy','sympy'])

        logger.debug('Done computing equations')
        
        #Compute the Jacobian (if required)
        if doJacobian:
            jac = sp.Matrix(dy).jacobian(yv).tolist()
            jacf = sp.lambdify([x,yv],jac,
                               modules=[{'M_P' : 1.22e19},labelsDict,'numpy','sympy'])
            logger.debug('Done computing Jacobian')
        else:
            jacf = None

        return rhsf,jacf

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
        S = self.normS*np.exp(r.y[-1,:])
        self.S = np.hstack((self.S,S))        
        #Store T-values
        NSvalues = r.y[-1,:]
        Tvalues = np.array([Tf(x,NSvalues[i],self.normS) for i,x in enumerate(r.t)])
        self.T = np.hstack((self.T,Tvalues))
        
        #Store the number and energy densities for each component:
        #(if the particle is coupled, use the equilibrium densities)
        for icomp,comp in enumerate(self.components):
            if not comp.active:
                n = np.array([np.nan]*len(r.t))
                rho = np.array([np.nan]*len(r.t))
            else:
                n = np.exp(r.y[icomp,:])*self.norm[icomp]
                rho = n*r.y[icomp+self.ncomp,:]
            comp.n = np.hstack((comp.n,n))
            comp.rho = np.hstack((comp.rho,rho))

    def printSummary(self,outFile=None):
        printSummary(self, outFile)
        
    def printData(self,outFile=None):
        printData(self, outFile)  

