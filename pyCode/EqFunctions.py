#!/usr/bin/env python3

"""

.. module:: AuxFuncs
    :synopsis: This module provides auxiliary functions 


:author: Andre Lessa <lessa.a.p@gmail.com>

"""

from pyCode.AuxFuncs import getFunctions
from sympy import Function
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

def nEQf(T,mass,dof):
    """Returns the equilibrium number density at temperature T. Returns zero for non-thermal components"""

    x = T/mass
    if x < 0.1:
        neq = mass**3*(x/(2*pi))**(3./2.)*exp(-1/x)*(1. + (15./8.)*x + (105./128.)*x**2)
    elif x < 1.5:
        neq = mass**3*x*besselk(2,1/x)/(2*pi**2)
    elif dof >= 0:
        neq = Zeta3*T**3/pi**2
    else:
        neq = (3./4.)*(Zeta3*T**3/pi**2)
        
    neq = neq*abs(dof)

    return neq

def dnEQdTf(T,mass,dof):
    """
    Returns the derivative of the equilibrium number density at temperature T
    with respect to T.
    """

    x = T/mass
    if x < 0.1:
        neq = mass**2*exp(-1/x)*(256. + 3*x*(288. + 5*x*(94. + 49*x)))/(512*sqrt(2.*x)*pi**(3./2.))
    elif x < 1.5:
        neq = mass**2*(besselk(1,1/x) + 3*x*besselk(2,1/x))/(2*x*pi**2)
    elif dof >= 0:
        neq = 3*Zeta3*T**2/pi**2
    else:
        neq = 3*(3./4.)*(Zeta3*T**2/pi**2)
        
    neq = neq*abs(dof)

    return neq

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

def drNeqdTf(T,massA,dofA,massB,dofB):
    """
    Returns the derivative of the ratio of equilibrium number densities
    with respect to T.
    If both species are Boltzmann suppressed (T << mass) and nearly degenerate,
    it provides a better approximation than taking the ratio of equilibrium densities,
    since the Boltzmann suppression exponentials partially cancel.
    """

    if massA == massB and dofA == dofB:
        return 1.

    x = T/massA
    y = T/massB
    if x < 0.1 and y < 0.1:
        r = exp(1/y-1/x)*(x-y)*sqrt(x*y)
        r *= (240*x*(-128+7*y*(y-16)) -128*(128 + 15*y*(16+7*y))+105*x**2*(-128 + y*(16+135*y)))
        if r:
            r = r/(T*x**3*(128+15*y*(16+7*y))**2)
    else:
        r = rNeqf(T,massA,dofA,massB,dofB) 
        if r:
            r *= dnEQdTf(T,massA,dofA)/nEQf(T,massA,dofA) - dnEQdTf(T,massB,dofB)/nEQf(T,massB,dofB)
            
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
        return (besselk(1,1/x)/besselk(2,1/x))*mass + 3.*x*mass
    else: #Non-relativistic
        return (1.-3.*x/2+15.*x**2/8)*mass + 3.*x*mass #Non-relativistic limit of bessel function ratio
    
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

class T(Function):
    """
    Defines a sympy function for the temperature.
    It will be undefined, unless evaluated numerically.
    Its derivative ignores changes in gSTARS, so it is defined to be T/3.
    """

    _imp_ = staticmethod(Tf)

    def fdiff(self, argindex=2):
        if argindex != 2:
            return 0
        else:
            return self/3

class nEQ(Function):
    """
    Defines a sympy function for the equilibrium number density
    as a function of temperature, mass and number of degrees of freedom.
    It will be undefined, unless evaluated numerically.
    """
    
    _imp_ = staticmethod(nEQf)

    def fdiff(self, argindex=1):
        if argindex == 1:
            return dnEQdT(self.args[0],self.args[1],self.args[2])
        else:
            return 0
    
class dnEQdT(Function):
    """
    Defines a sympy function for the derivative of the
     equilibrium number density with respect to temperature.
    It will be undefined, unless evaluated numerically.
    """    
    
    _imp_ = staticmethod(dnEQdTf)
    
class Pn(Function):
    """
    Defines a sympy function for the ratio of pressure and number density
    for a component given the ratio of its energy and number densities and
    the mass.
    It will be undefined, unless evaluated numerically.
    """
    
    _imp_ = staticmethod(Pnf)

    def fdiff(self, argindex=1):
        if argindex == 1:
            return dPndR(self.args[0],self.args[1])
        else:
            return 0
    
class dPndR(Function):
    """
    Defines a sympy function for the derivative of the ratio of pressure and number density
    with respect to the  ratio of the energy and number densities.
    It will be undefined, unless evaluated numerically.
    """    
    
    _imp_ = staticmethod(dPndRf)    

class rNeq(Function):
    """
    Defines a sympy function for the ratio of equilibrium number densities
    as a function of temperature, mass and number of degrees of freedom
    It will be undefined, unless evaluated numerically.
    """
    
    _imp_ = staticmethod(rNeqf)

    def fdiff(self, argindex=1):
        if argindex == 1:
            return drNeqdT(self.args[0],self.args[1],self.args[2],
                            self.args[3],self.args[4])
        else:
            return 0
    
class drNeqdT(Function):
    """
    Defines a sympy function for the derivative of the
    ratio of equilibrium number densities with respect to temperature.
    It will be undefined, unless evaluated numerically.
    """    
    
    _imp_ = staticmethod(drNeqdTf)

class rEQ(Function):
    """
    Defines a sympy function for the ratio of equilibrium energy 
    and number densities at temperature T,
    assuming chemical potential = 0.
    It will be undefined, unless evaluated numerically.
    """
    
    _imp_ = staticmethod(rEQf)

    def fdiff(self, argindex=1):
        if argindex == 1:
            return drEQdT(self.args[0],self.args[1],self.args[2])
        else:
            return 0
    
class drEQdT(Function):
    """
    Defines a sympy function for the derivative of the
    ratio of equilibrium energy and number densities 
    with respect to the temperature T. Assumes zero chemical potential.
    It will be undefined, unless evaluated numerically.
    """    
    
    _imp_ = staticmethod(drEQdTf)


