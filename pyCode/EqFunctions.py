#!/usr/bin/env python3

"""

.. module:: AuxFuncs
    :synopsis: This module provides auxiliary functions 


:author: Andre Lessa <lessa.a.p@gmail.com>

"""

from pyCode.AuxFuncs import getFunctions
from sympy import Function,exp,sqrt,log, pi
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
import warnings
Tmin,Tmax = 1e-15,1e5 #min and max values for evaluating gSTAR

warnings.filterwarnings('error')

Tfunc,gSTARf,gSTARSf = getFunctions()
        
class gSTAR(Function):
    """
    Defines a sympy function for the effective number of relativistic
    degrees of freedom for the energy density.
    It will be undefined, unless evaluated numerically.
    Its derivative is always taken to be zero, since it is a
    smooth function.
    """

    _imp_ = staticmethod(gSTARf)

    def fdiff(self, argindex=1):
        return 0

class gSTARS(Function):
    """
    Defines a sympy function for the effective number of relativistic
    degrees of freedom for the entropy.
    It will be undefined, unless evaluated numerically.
    Its derivative is always taken to be zero, since it is a
    smooth function.
    """

    _imp_ = staticmethod(gSTARSf)

    def fdiff(self, argindex=1):
        return 0

def H(T, rhoTot):
    """Compute the Hubble parameter, given the temperature and total energy density"""
    
    MP = 1.22*10**19
    
    rhoRad = (pi**2/30)*gSTAR(T)*T**4  # thermal bath's energy density    
    rho = rhoRad+rhoTot
    H = sqrt(8*pi*rho/3)/MP
    
    return H

def getTemperature(x,NS,S0=1.):
    
    xmin = log((2*pi**2/45.)*gSTARS(Tmin)*Tmin**3)
    xmax = log((2*pi**2/45.)*gSTARS(Tmax)*Tmax**3)
    xeff = NS + log(S0) - 3.*x
        
    if xeff < xmin:  #For T < Tmin, g* is constant
        return ((45./(2*pi**2))*exp(xeff)/gSTARS(Tmin))**(1./3.)
    elif xeff > xmax: #For T > Tmax, g* is constant
        return ((45./(2*pi**2))*exp(xeff)/gSTARS(Tmax))**(1./3.)
    else:    
        return float(Tfunc(xeff))

class T(Function):
    """
    Defines a sympy function for the temperature.
    It will be undefined, unless evaluated numerically.
    Its derivative ignores changes in gSTARS, so it is defined to be T/3.
    """

    _imp_ = staticmethod(getTemperature)

    def fdiff(self, argindex=2):
        if argindex != 2:
            return 0
        else:
            return self/3

