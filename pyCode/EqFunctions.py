#!/usr/bin/env python3

"""

.. module:: AuxFuncs
    :synopsis: This module provides auxiliary functions 


:author: Andre Lessa <lessa.a.p@gmail.com>

"""

from pyCode.AuxFuncs import getFunctions
from numpy import pi,sqrt,exp,log
from numpy.polynomial.polynomial import polyval,polyder
from scipy.special import kn as besselk
from scipy.special import zetac
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
import warnings
Tmin,Tmax = 1e-15,1e5 #min and max values for evaluating gSTAR
Zeta3 = zetac(3.) + 1.

warnings.filterwarnings('error')

Tfunc,gSTARf,gSTARSf = getFunctions()
#Dummy definitions:
gSTAR = gSTARf
gSTARS = gSTARSf

def Tf(x,NS,S0=1.):
    
    xmin = log((2*pi**2/45.)*gSTARSf(Tmin)*Tmin**3)
    xmax = log((2*pi**2/45.)*gSTARSf(Tmax)*Tmax**3)
    xeff = NS + log(S0) - 3.*x
        
    if xeff < xmin:  #For T < Tmin, g* is constant
        return ((45./(2*pi**2))*exp(xeff)/gSTARSf(Tmin))**(1./3.)
    elif xeff > xmax: #For T > Tmax, g* is constant
        return ((45./(2*pi**2))*exp(xeff)/gSTARSf(Tmax))**(1./3.)
    else:    
        return float(Tfunc(xeff))

#Define relevant derivatives:    
# dTfdx = lambda x,NS,S0: derivative(Tf,x,args=(NS,S0),dx=1e-5)
#Approximation (neglects variation in gSTARS):
def dTfdx(x,NS,S0):
    """
    Returns the derivative of the temperature T
    with respect to x, assuming that the variation of gSTARS is negligible.
    """   

    return -Tf(x,NS,S0)

def nEQf(T,mass,dof):
    """
    Returns the equilibrium number density at temperature T. Returns zero for non-thermal components
    """

    a = 15./8.
    b = 105./129.
    if mass:
        x = T/mass
    else:
        x = 1e3 #Use ultra relativistic limit for massless particles

    if x < 0.1:
        neq = mass**3*(x/(2*pi))**(3./2.)*exp(-1/x)*(1. + a*x + b*x**2) #Expansion at low x
    elif x < 1.5:
        neq = mass**3*x*besselk(2,1/x)/(2*pi**2)
    elif dof >= 0:
        neq = Zeta3*T**3/pi**2
    else:
        neq = (3./4.)*(Zeta3*T**3/pi**2)
        
    neq = neq*abs(dof)

    return neq

def dnEQfdT(T,mass,dof):
    """
    Returns the derivative of the equilibrium number density at temperature T
    with respect to T.
    """

    a = 15./8.
    b = 105./129.
    if mass:
        x = T/mass
    else:
        x = 1e3 #Use ultra relativistic limit for massless particles    
    
    if not x:
        return 0.   
    if x < 0.1:
        dneq = mass**2*exp(-1/x)*(2 + (3+2*a)*x + (5*a+2*b)*x**2 + 7*b*x**3)
        dneq = dneq/(2*sqrt(x)*(2*pi)**(3/2))
    elif x < 1.5:
        dneq = mass**2*(besselk(1,1/x) + 3*x*besselk(2,1/x))/(2*x*pi**2)
    elif dof >= 0:
        dneq = 3*Zeta3*T**2/pi**2
    else:
        dneq = 3*(3./4.)*(Zeta3*T**2/pi**2)
        
    dneq = dneq*abs(dof)

    return dneq

def dLnEQfdT(T,mass):
    """
    Returns the derivative of the equilibrium number density at temperature T
    with respect to T divided by the equilibrium number density.
    """

    a = 15./8.
    b = 105./129.
    if mass:
        x = T/mass
    else:
        x = 1e3 #Use ultra relativistic limit for massless particles    

    if not x or not T:
        return 0.    
    if x < 0.1:
        dLneq = 1/x**2 + 3/(2*x) + (a+2*b*x)/(1+x*(a+b*x))
        dLneq = dLneq/mass
    elif x < 1.5:
        dLneq = 3/x + (1/x**2)*besselk(1,1/x)/besselk(2,1/x)
        dLneq = dLneq/mass
    else:
        dLneq = 3/T

    return dLneq

def Pnf(R,mass):
    """
    Computes the ratio of pressure and number density
    for a component given the ratio of its energy and number densities and
    the mass.
    """

    r = R/mass
    Pnrel = (R/3)  # Relativistic pressure

    if r > 11.4:
        return Pnrel  # Ultra relativistic limit
    if r <= 1:
        return 0.  # Ultra non-relativistic limit
    
    # Expansion coefficients for relativistic/non-relativistic transition    
    aV = [0.,2./3.,-0.345998, 0.234319, -0.0953434, 
          0.023657, -0.00360707, 0.000329645, 
          -0.0000165549, 3.51085e-7]
    x = r-1 #expansion in terms of r-1
    Pnnonrel = mass*polyval(x,aV)
        
    return Pnnonrel

def dPndRf(R,mass):
    """
    Computes the derivative of the ratio of pressure and number density
    with respect to the  ratio of the energy and number densities.
    """
    
    r = R/mass
    if r > 11.4:
        return 1./3.
    if r <= 1.:
        return 0.
    
    aV = [0.,2./3.,-0.345998, 0.234319, -0.0953434, 
          0.023657, -0.00360707, 0.000329645, 
          -0.0000165549, 3.51085e-7]     
    daV = polyder(aV) #coefficients of the derivative
    x = r-1
    dPnnonrel = polyval(x,daV)

    return dPnnonrel

def rNeqf(T,massA,dofA,massB,dofB):
    """
    Returns the ratio of equilibrium number densities at temperature T (nEQA/nEQB).
    If both species are Boltzmann suppressed (T << mass) and nearly degenerate,
    it provides a better approximation than taking the ratio of equilibrium densities,
    since the Boltzmann suppression exponentials partially cancel.
    """

    if massA == massB and dofA == dofB:
        return 1.

    x = T/massA
    y = T/massB
    if x < 0.1 and y < 0.1:
        r = (y/x)**(3./2.)*exp(1/y-1/x)
        r *= (1. + (15./8.)*x + (105./128.)*x**2)/(1. + (15./8.)*y + (105./128.)*y**2)
        r *= abs(dofA)/abs(dofB)
    else:
        r = nEQf(T,massA,dofA)
        if r:
            r = r/nEQf(T,massB,dofB)
            
    return r

def rEQf(T,mass,dof):
    """
    Returns the ratio of equilibrium energy and number densities at temperature T,
    assuming chemical potential = 0.
    """

    x = T/mass
    
    if x > 1.675:   #Ultra relativistic
        if dof < 0:
            return (7./6.)*pi**4*T/(30.*Zeta3)  #Fermions
        else:
            return pi**4*T/(30.*Zeta3)    #Bosons
    elif x > 1e-2:  #Non-relativistic/relativistic transition
        return mass*(besselk(1,1/x)/besselk(2,1/x) + 3.*x)
    else: #Non-relativistic
        return mass*((1.-3.*x/2+15.*x**2/8) + 3.*x) #Non-relativistic limit of bessel function ratio
    
def drEQdTf(T,mass,dof):
    """
    Returns the derivative of the ratio of equilibrium energy and number densities 
    with respect to the temperature T. Assumes zero chemical potential.
    """

    x = T/mass
    
    if x > 1.675:   #Ultra relativistic
        return rEQf(T,mass,dof)/T
    elif x > 1e-2:  #Non-relativistic/relativistic transition
        r = 1./(2*x**2*besselk(2,1/x)**2)
        r *= (-besselk(2,1/x)*(besselk(0,1/x)+(1+6*x**2)*besselk(2,1/x))
              + besselk(1,1/x)*(besselk(1,1/x)+besselk(3,1/x)))
        return r
    else: #Non-relativistic
        return (3./4.)*(2+5*x)    
