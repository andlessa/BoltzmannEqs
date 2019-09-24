#!/usr/bin/env python3

"""

.. module:: AuxFuncs
    :synopsis: This module provides auxiliary functions 


:author: Andre Lessa <lessa.a.p@gmail.com>

"""

import os
from scipy import integrate, interpolate
from sympy import Function
from numpy import logspace, log10,exp,sqrt,log, pi 
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
import warnings
import pickle
Tmin,Tmax = 1e-15,1e5 #min and max values for evaluating gSTAR

warnings.filterwarnings('error')


class interp1d_picklable:
    """
    class wrapper for piecewise linear function. Required for pickling a interp1d result.
    """
    def __init__(self, xi, yi, **kwargs):
        self.xi = xi
        self.yi = yi
        self.args = kwargs
        self.f = interpolate.interp1d(xi, yi, **kwargs)

    def __call__(self, xnew):
        return self.f(xnew)

    def __getstate__(self):
        return self.xi, self.yi, self.args

    def __setstate__(self, state):
        self.f = interpolate.interp1d(state[0], state[1], **state[2])
        
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

def Hfunc(T, rhoTot):
    """Compute the Hubble parameter, given the temperature and total energy density"""
    
    MP = 1.22*10**19
    
    rhoRad = (pi**2/30)*gSTAR(T)*T**4  # thermal bath's energy density    
    rho = rhoRad+rhoTot
    H = sqrt(8*pi*rho/3)/MP
    
    return H

def getFunctions(pclFile):
    """
    Computes the g*(T), g*s(T) and temperature functions and saves
    them to a pickle file. Ignores all BSM effects to these functions
    :param pclFile: Name of pickle file to dump the functions
    """
    
    logger.info("Computing auxiliary functions. This calculation is done only once and the results will be stored in %s.\n" %pclFile)
    
    #Get points to evaluate gSTAR
    Tpts = logspace(log10(Tmin),log10(Tmax),num=2000,endpoint=False)
    #Evaluate gSTAR and gSTARS at these points
    gSTARpts = [gSTARexact(T) for T in Tpts]
    gSTARSpts = [gSTARSexact(T) for T in Tpts]
    #Get interpolating functions:
    gSTAR = interp1d_picklable(Tpts,gSTARpts,fill_value = (gSTARpts[0],gSTARpts[-1]),
                               bounds_error=False)
    gSTARS = interp1d_picklable(Tpts,gSTARSpts,fill_value = (gSTARSpts[0],gSTARSpts[-1]),
                               bounds_error=False)
    #Evaluate (2*pi^2/45)*gstarS(T)*T^3 at these points:
    fpts = [log((2*pi**2/45.)*gSTARS(T)*T**3) for T in Tpts]
    #Get inverse function to compute temperature from 
    Tfunc =  interp1d_picklable(fpts,Tpts,fill_value='extrapolate')    
    f = open(pclFile,'wb')
    pickle.dump(gSTAR,f)
    pickle.dump(gSTARS,f)
    pickle.dump(Tfunc,f)
    f.close()

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

def gstarFunc(x, dof):
    """
    Auxiliary function to compute the contribution from a single particle to gSTAR.\
    x = mass/T, dof = number of degrees of freedom (positive/negative for bosons/fermions)
    """
            
    if x > 20.: res = 0.  # Particle has decoupled
    elif x > 10.**(-2):  # Near decoupling
        ep = -float(dof / abs(dof))
        epsilon = 0.01  # To avoid limits on end points
        a = 0. + epsilon
        b = (1. / x) * (1. - epsilon / 100.)
        res = integrate.romberg(lambda y: sqrt(1. - y ** 2 * x ** 2) / (y ** 5 * (exp(1. / y) + ep)), a, b,
                                tol=0., rtol=0.01, divmax=100, vec_func=False)
    else:
        if dof < 0: res = 5.6822  # Fully relativistic/coupled
        elif dof > 0: res = 6.49394
                        
    return res*abs(dof)*0.15399  # Result

def gSTARexact(T, interpol=True):
    """
    Computes exactly the number of relativistic degrees of freedom in thermal equilibrium with the thermal
    bath.
    interpol turns on/off the interpolation around the QCD phase trasition region.
    """

    gstar = 0.
# Define SM masses and degrees of freedom
    MassesGauge = {"W" : 80., "Z" : 91., "A" : 0.}
    DoFGauge = {"W" : 6, "Z" : 3, "A" : 2}
    MassesLeptons = {"electron" : 0.51 * 10.**(-3), "muon" : 0.1056, "tau" : 1.77, "neutrino" : 0.}
    DoFLeptons = {"electron" :-4, "muon" :-4, "tau" :-4, "neutrino" :-6}
    if T < 0.25:  # After QCD phase transition        
        MassesHadrons = {"pion" : 0.14, "eta" : 0.55, "rho" : 0.77, "omega" : 0.78, "kaon" : 0.5}
        DoFHadrons = {"pion" : 4, "eta" : 2, "rho" : 6, "omega" : 6, "kaon" : 4}
    else:  # Before QCD phase transition
        MassesHadrons = {"u" : 3.*10 ** (-3), "d" : 5.*10 ** (-3), "s" : 0.1, "c" : 1.3, "b" : 4.2, "t" : 173.3,
                         "g" : 0.}
        DoFHadrons = {"u" :-12, "d" :-12, "s" :-12, "c" :-12, "b" :-12, "t" :-12, "g" : 16}
    MassesSM = dict(list(MassesGauge.items()) + list(MassesLeptons.items()) + list(MassesHadrons.items()))
    DoFSM = dict(list(DoFGauge.items()) + list(DoFLeptons.items()) + list(DoFHadrons.items()))
# Add up SM degrees of freedom     
    for part in MassesSM:
        gstar += gstarFunc(MassesSM[part] / T, DoFSM[part])

# Correct for neutrino decoupling:
    if T <= 5.*10.**(-4):
        gstar += (-1. + (4. / 11.) ** (4. / 3.)) * gstarFunc(MassesSM["neutrino"] / T, DoFSM["neutrino"]) 
     
# Smooth discontinuous transitions:
    if interpol:
# QCD phase transition:
        finter = None
        if 0.15 < T < 0.3 and interpol:
            Tpts = [0.15, 0.3]
            gpts = [gSTARexact(Tpt, False) for Tpt in Tpts]        
            finter = interpolate.interp1d(Tpts, gpts, kind='linear')            
# Neutrino decoupling            
        elif  2.*10.**(-4) < T < 6.*10.**(-4):
            Tpts = [2.*10.**(-4), 6.*10 ** (-4)]
            gpts = [gSTARexact(Tpt, False) for Tpt in Tpts]        
            finter = interpolate.interp1d(Tpts, gpts, kind='linear')
# Replace gstar value by interpolation:
        if finter:
            gstar = finter(T)
    
    return gstar

def gSTARSexact(T):
    """
    Computes the number of relativistic degrees of freedom for computing the entropy density,\
    including the full MSSM spectrum, except for the lightest neutralino.
    """
 
    if T >= 10.**(-3):
        return gSTARexact(T)
    else:  # Correct for neutrino decoupling:        
        return gSTARexact(T) - (7. / 8.) * 6.*(4. / 11.) ** (4. / 3.) + (7. / 8.) * 6.*(4. / 11.)


#Load auxiliary (pre-computed) functions:
if not os.path.isfile('gFunctions.pcl'):
    getFunctions('gFunctions.pcl')
 
f = open('gFunctions.pcl','rb')
logger.info("Loading aux functions. Ignoring BSM corrections to g* and g*_S")
gSTARf = pickle.load(f)
gSTARSf = pickle.load(f)
Tfunc = pickle.load(f)
f.close()
