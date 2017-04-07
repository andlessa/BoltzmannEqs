#!/usr/bin/env python

"""

.. module:: COsc
    :synopsis: This module defines all the oscillation amplitudes for the components at temperature T


:synopsis: This module defines all the oscillation amplitudes for the components at temperature T
:author: Andre Lessa <lessa.a.p@gmail.com>

"""


IDs = ['axino', 'saxion', 'saxionCO', 'axion', 'axionCO', 'gravitino', 'neutralino']

import modelParameters
from parameters import Pi
from math import log, exp
import Mass
from AuxFuncs import gSTARS

def axino(T):
    """Axino coherent oscillation amplitude as a function of temperature """
    
    return 0.

def saxion(T):
    """Saxion coherent oscillation amplitude as a function of temperature """
    
    return 0.

def saxionCO(T):
    """Saxion coherent oscillation amplitude as a function of temperature """
    
    TR = modelParameters.TR
    sI = modelParameters.sI
    
    r_saxion1 = 1.9*10.**(-8)*(2.*Pi**2*gSTARS(TR)*TR**3/45.)*(TR/10.**5)*(sI/10.**12)**2
    r_saxion2 = Mass.saxion(T)**2*sI**2/2. 
    
    return min(r_saxion1,r_saxion2)    
    

def axion(T):
    """Axion coherent oscillation amplitude as a function of temperature """
    
    return 0.
         
def axionCO(T):
    """Axion coherent oscillation amplitude as a function of temperature """
    
    fa = modelParameters.fa
    thetaI = modelParameters.thetaI
    
    if modelParameters.modelType == 'DFSZ': fa = fa/6.
    Ftheta = (log(exp(1.)/(1.-thetaI**2/Pi**2)))**(7./6.)
    r_axion = 1.44*fa**2*Mass.axion(T)**2*thetaI**2*Ftheta/2.        
    return r_axion

def gravitino(T):
    """Gravitino coherent oscillation amplitude as a function of temperature """
    
    return 0.

def neutralino(T):
    """Neutralino1 coherent oscillation amplitude as a function of temperature """
    
    return 0.  

