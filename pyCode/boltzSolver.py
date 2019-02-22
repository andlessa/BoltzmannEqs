#!/usr/bin/env python

"""

.. module:: boltzSolver
    :synopsis: This module contains the main methods for solving the Boltzmann equations 

:synopsis: This module contains the main methods for solving the Boltzmann equations
:author: Andre Lessa <lessa.a.p@gmail.com>

"""

from pyCode.boltzEqs import BoltzEqs
from pyCode.AuxFuncs import gSTARS, getTemperature
from math import log,exp,pi
import logging
import random, time
import numpy as np
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)
random.seed('myseed')


def Evolve(compList,T0,TF,omegaErr=0.01):
    """Evolve the components in component list from the re-heat temperature T0 to TF
    For simplicity we set  R0 = s0 = 1 (with respect to the notes).
    Returns a list with components, where the evolveVars hold the evolution of each component.
    The last element of the list is a simple array with the evolution of NS.
    omegaErr is the approximate relative error for Omega h^2.    
    """

#Sanity checks
    if not goodCompList(compList,T0): return False
    
    t0 = time.time()
#Compute initial conditions    
    for comp in compList:
        comp.setInitialCond(T0)
    x0 = 0.      # Initial condition for log(R/R0)
    y0 = [comp.evolveVars["N"] for comp in compList]  #Initial conditions for log(n/s0)
    y0 += [comp.evolveVars["R"] for comp in compList] #Initial conditions for log(rho/n)
    S = (2*pi**2/45)*gSTARS(T0)*T0**3
    y0.append(log(S))  #Initial condition for log(S/S0)
    sw = [comp.active for comp in compList]
    logger.info("Initial conditions computed in %s s" %(time.time()-t0))
    t0 = time.time()
    
#Solve differential equations:
#First call with large errors to estimate the size of the solution
    xf = 50.
#     rtol = omegaErr
#     atol = [omegaErr]*len(y0)    
    boltz_eqs = BoltzEqs(compList,x0,y0,sw) #Define equations and set initial conditions
    
    boltz_eqs.solver.set_initial_value(y0,x0)
    boltz_eqs.solver.set_integrator('lsoda',method='bdf')    
    x_vals = np.linspace(x0, xf, 100)
    vals = []
    for xv in x_vals[1:]:
        x_init,y_init = boltz_eqs.solver.t,boltz_eqs.solver.y
        nsteps = 10
        boltz_eqs.solver.set_integrator('lsoda',nsteps=nsteps)
        boltz_eqs.solver.set_initial_value(y_init,x_init)        
        print('x0=',x_init,'xf=',xv,'nsteps=',nsteps,boltz_eqs.solver.get_return_code())
        while (boltz_eqs.solver.get_return_code() != 2) and nsteps < 10000:
            try:
                boltz_eqs.solver.set_integrator('lsoda',method='bdf',nsteps=nsteps)
                boltz_eqs.solver.set_initial_value(y_init,x_init)
                boltz_eqs.solver.integrate(xv)
                vals.append(np.array([boltz_eqs.solver.t] + boltz_eqs.solver.y.tolist()))                
            except Exception as e:
                nsteps *= 5
                print(e)
                print(nsteps)
        if boltz_eqs.solver.get_return_code() != 2:
            break
        

    return np.array(vals)
    
    
#     while boltz_eqs.solver.successful() and boltz_eqs.solver.t < xf:
#         print('x=',boltz_eqs.solver.t)
#         boltz_eqs.solver.integrate(boltz_eqs.solver.t + x_step)
#         xvals.append(boltz_eqs.solver.t)
#         yvals.append(boltz_eqs.solver.y)    
    
#     y = boltz_eqs.solver.integrate(xf)
#     logger.info("First pass at solving Boltzmann equations done in %s s" %(time.time()-t0))
#     t0 = time.time()

#     for comp in compList: comp.evolveVars = {'T' : [], 'R' : [], 'rho' : [], 'n' : []}
#     for ipt,ypt in enumerate(y):
#         NS = ypt[-1]
#         T = getTemperature(x[ipt],NS)
#         if T < max(10.**(-7),TF): continue     #Do not keep points above TF or after matter domination
#         for icomp,comp in enumerate(compList):        
#             comp.evolveVars['T'].append(T)
#             comp.evolveVars['R'].append(exp(x[ipt]))
#             n = exp(ypt[icomp])
#             if comp.Type == 'CO': rho = n*comp.mass(T)
#             else: rho = n*ypt[icomp + len(compList)]
#             if T < comp.Tdecay or (comp.Type == 'CO' and T > comp.Tosc): n = rho = 0.            
#             comp.evolveVars['rho'].append(rho)
#             comp.evolveVars['n'].append(n)
#     return True
#     return np.array(xvals),np.array(yvals)

        
                
def goodCompList(compList,T0):
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
    
    
