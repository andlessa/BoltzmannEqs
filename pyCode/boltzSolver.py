#!/usr/bin/env python

"""

.. module:: boltzSolver
    :synopsis: This module contains the main methods for solving the Boltzmann equations 

:synopsis: This module contains the main methods for solving the Boltzmann equations
:author: Andre Lessa <lessa.a.p@gmail.com>

"""

from pyCode.boltzEqs import BoltzEqs
from pyCode.AuxFuncs import gSTARS, getTemperature
from numpy import log,exp,pi
import numpy as np
from scipy import integrate
import logging
import random, time
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
    tvals = np.linspace(x0,xf,100)
    boltz_eqs = BoltzEqs(compList,x0,y0,sw) #Define equations and set initial conditions
    r = integrate.solve_ivp(boltz_eqs.rhs,t_span=(x0,xf),y0=y0,
                            t_eval=tvals,method='BDF',
                            events=boltz_eqs.events)
    xvals = r.t
    yvals = r.y
    while r.t[-1] < xf and r.status >= 0:
        if r.t_events:
            x0,y0,sw = boltz_eqs.handle_events(r)
            boltz_eqs.updateValues(x0, y0, sw)
            tvals = [t for t in tvals[:] if t > x0]
            r = integrate.solve_ivp(boltz_eqs.rhs,t_span=(x0,xf),y0=y0,
                                    t_eval=tvals,method='BDF',
                                    events=boltz_eqs.events)
            xvals = np.concatenate((xvals,r.t))
            yvals = np.concatenate((yvals,r.y),axis=1)            
    
    if r.status < 0:
        logger.error(r.message)

    return xvals,yvals
    
        
                
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
    
    
