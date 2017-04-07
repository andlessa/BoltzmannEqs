#!/usr/bin/env python

"""

.. module:: boltzSolver
    :synopsis: This module contains the main methods for solving the Boltzmann equations 

:synopsis: This module contains the main methods for solving the Boltzmann equations
:author: Andre Lessa <lessa.a.p@gmail.com>

"""

from boltzEqs import BoltzEqs
from assimulo.solvers import CVode
from parameters import Pi
from AuxFuncs import gSTARS, getTemperature
from math import log,exp
import logging
import random
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
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
            
#Compute initial conditions    
    for comp in compList: comp.setInitialCond(T0)
    x0 = 0.      # Initial condition for log(R/R0)
    y0 = [comp.evolveVars["N"] for comp in compList]  #Initial conditions for log(n/s0)
    y0 += [comp.evolveVars["R"] for comp in compList] #Initial conditions for log(rho/n)
    S = (2*Pi**2/45)*gSTARS(T0)*T0**3
    y0.append(log(S))  #Initial condition for log(S/S0)
    sw = [comp.active for comp in compList]
    
#Solve differential equations:
#First call with large errors to estimate the size of the solution
    xf = 50.
    rtol = omegaErr
    atol = [omegaErr]*len(y0)    
    boltz_eqs = BoltzEqs(compList,x0,y0,sw) #Define equations and set initial conditions
    y = mySolve(xf,boltz_eqs,rtol,atol)[1]
#Second call with proper relative/absolute errors
    maxy = max([abs(yy) for yy in y[-1][:len(compList)]])
    rtol = omegaErr/(2.*maxy)  #Re-scale error to obtain relic densities with precision omegaErr
    atol = [omegaErr/2.]*len(compList) + [omegaErr*abs(yy) for yy in y[-1][len(compList):]]
    atol = [max(xx,0.005) for xx in atol]
    boltz_eqs = BoltzEqs(compList,x0,y0,sw) #Define equations and set initial conditions
    x,y = mySolve(xf,boltz_eqs,rtol,atol,verbosity=50)

#Store the solutions:    
    for comp in compList: comp.evolveVars = {'T' : [], 'R' : [], 'rho' : [], 'n' : []}
    for ipt,ypt in enumerate(y):
        NS = ypt[-1]
        T = getTemperature(x[ipt],NS)
        if T < max(10.**(-7),TF): continue     #Do not keep points above TF or after matter domination
        for icomp,comp in enumerate(compList):        
            comp.evolveVars['T'].append(T)
            comp.evolveVars['R'].append(exp(x[ipt]))
            n = exp(ypt[icomp])
            if comp.Type == 'CO': rho = n*comp.mass(T)
            else: rho = n*ypt[icomp + len(compList)]
            if T < comp.Tdecay or (comp.Type == 'CO' and T > comp.Tosc): n = rho = 0.            
            comp.evolveVars['rho'].append(rho)
            comp.evolveVars['n'].append(n)

    return True

        
def mySolve(xf,boltz_eqs,rtol,atol,verbosity=50):
    """Sets the main options for the ODE solver and solve the equations. Returns the
    array of x,y points for all components.
    If numerical instabilities are found, re-do the problematic part of the evolution with smaller steps"""
        
    boltz_solver = CVode(boltz_eqs)  #Define solver method
    boltz_solver.rtol = rtol
    boltz_solver.atol = atol
    boltz_solver.verbosity = verbosity
    boltz_solver.linear_solver = 'SPGMR'
    boltz_solver.maxh = xf/300.
    xfinal = xf
    xres = []
    yres = []
    sw = boltz_solver.sw[:]
    while xfinal <= xf:
        try:
            boltz_solver.re_init(boltz_eqs.t0,boltz_eqs.y0)
            boltz_solver.sw = sw[:]
            x,y = boltz_solver.simulate(xfinal)
            xres += x
            for ypt in y: yres.append(ypt)
            if xfinal == xf: break   #Evolution has been performed until xf -> exit            
        except Exception,e:
            print e
            if not e.t or 'first call' in e.msg[e.value]:
                logger.error("Error solving equations:\n "+str(e))
                return False
            xfinal = max(e.t*random.uniform(0.85,0.95),boltz_eqs.t0+boltz_solver.maxh)  #Try again, but now only until the error
            logger.warning("Numerical instability found. Restarting evolution from x = "
                           +str(boltz_eqs.t0)+" to x = "+str(xfinal))
            continue
        xfinal = xf  #In the next step try to evolve from xfinal -> xf
        sw = boltz_solver.sw[:]
        x0 = float(x[-1])
        y0 = [float(yval) for yval in y[-1]]
        boltz_eqs.updateValues(x0,y0,sw)

    
    return xres,yres
                            
                
def goodCompList(compList,T0):
    """Checks if the list of components satisfies the minimum requirements"""
    
    if type(compList) != type(list()):
        logger.error("Input must be a list of Component objects")
        return False
    for comp1 in compList:
        BRs1 = comp1.getBRs(T0)
        for comp2 in compList:            
            if comp1 == comp2: continue
            if not comp2.ID in BRs1.getAllFinalStates(): continue
            if comp1.width(T0) < comp2.width(T0):
                logger.error("Improper decay from %s to %s" %(comp1.ID,comp2.ID))
                return False

    return True
    
    
