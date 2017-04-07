#!/usr/bin/env python

"""

.. module:: Mass
    :synopsis: This module defines all the mass functions for the relevant components\
    It imports the model parameters from modelParameters.\
    The function names must match the particle IDs
    
:synopsis: This module defines all the mass functions for the relevant components\
    It imports the model parameters from modelParameters.\
    The function names must match the particle IDs    
:author: Andre Lessa <lessa.a.p@gmail.com>

"""


IDs = ['axino', 'saxion', 'saxionCO', 'axion', 'axionCO', 'gravitino', 'neutralino']

import modelParameters
from parameters import LQCD

def axino(T):
    """Axino mass as a function of temperature """
    
    return abs(modelParameters.mass_axino)

def saxion(T):
    """Saxion mass as a function of temperature """
    
    return abs(modelParameters.mass_saxion)

def saxionCO(T):
    """Saxion mass as a function of temperature """
    
    return abs(saxion(T))    
    

def axion(T):
    """Axion mass as a function of temperature """

    fa = modelParameters.fa
    if modelParameters.modelType == 'DFSZ': fa = fa/6.
            
    return (6.2*10**(-3))*min(1.,0.018*(LQCD/T)**4)/fa
         
def axionCO(T):
    """Axion mass as a function of temperature """
    
    return axion(T)    

def gravitino(T):
    """Gravitino mass as a function of temperature """
    
    return abs(modelParameters.mass_gravitino)

def neutralino(T):
    """Neutralino1 mass as a function of temperature """
    
    return abs(modelParameters.Masses[1000022]) 

