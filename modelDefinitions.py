#!/usr/bin/env python

"""

.. module:: modelParameters
    :synopsis: This module defines all relevant input parameters for the model and some methods to obtain \
    other relevant information
    
:synopsis: This module defines all relevant input parameters for the model and some methods to obtain \
    other relevant information (spectrum, decays, sigma.v, etc)    
:author: Andre Lessa <lessa.a.p@gmail.com>

"""

from pyCode.AuxDecays import DecayList, Decay
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


mDM = 100.
mMediator = 500.
yCoupling = (1.5e-13)
alpha = 1./137.


def MediatorDecays(T):
    """
    Defines the mediator decays.
    
    :param T: temperature
    :return: DecayList object with the mediator decays and width
    """
    
    decays = DecayList()
    decayToDM = Decay(instate='Mediator',fstates=['DM','radiation'],br=1.)
    decays.addDecay(decayToDM)
    decays.Xfraction = 0.5 #
    decays.width = yCoupling**2*mMediator/(8.*3.14)
    
    return decays

def MediatorSigmaV(T):
    """
    Defines the mediator annihilation cross-section as a function of temperature.
    
    :param T: temperature
    :return: thermally averaged cross-section
    """
    
    return alpha/mMediator**2

def DMSigmaV(T):
    """
    Defines the dark matter annihilation cross-section as a function of temperature.
    
    :param T: temperature
    :return: thermally averaged cross-section
    """
    
    return yCoupling**4/mDM**2

