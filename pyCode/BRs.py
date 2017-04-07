#!/usr/bin/env python

"""

.. module:: BRs
    :synopsis: This module defines all the BRs and widths for the relevant components\
    It imports the model parameters from modelParameters.\
    The function names must match the particle IDs

:synopsis: This module defines all the BRs and widths for the relevant components\
    It imports the model parameters from modelParameters.\
    The function names must match the particle IDs
:author: Andre Lessa <lessa.a.p@gmail.com>

"""


import axinoDecays, saxionDecays, n1Decays, gravitinoDecays
import modelParameters
from AuxFuncs import memoize
from AuxDecays import DecayList

IDs = ['axino', 'saxion', 'saxionCO', 'axion', 'axionCO', 'gravitino', 'neutralino', 'radiation']
susyQ = modelParameters.slhadata['MSOFT'].q  #SUSY scale


def compressBRs(In_BRs,T):
    """Checks if the final states (B) have a lifetime shorter than the parent (A). If they do, compress the decay,
    so there are no improper decays. Explicitly, if A -> B(->D+E) + C  and\
    lifetime(B) < lifetime(A), then BR(A->D + E + C) = BR(A-> B + C)* BR(B-> D + E)."""
    
    if In_BRs.width == 0.: return In_BRs #Do nothing for stable particles
    
    Out_BRs = In_BRs
    earlyDecays = {'radiation': DecayList()}
    for ID in In_BRs.getAllFinalStates():       
        if not ID in earlyDecays:
            earlyDecays[ID] = eval(ID+'BRs(T)')
        if earlyDecays[ID].width <= In_BRs.width:  earlyDecays.pop(ID)    

    unstableIDs = set(Out_BRs.getAllFinalStates()).intersection(set(earlyDecays.keys()))
    while unstableIDs: 
        newDecayList = DecayList()
        newDecayList.width = Out_BRs.width
        newDecayList.Xfraction = Out_BRs.Xfraction
        ID = unstableIDs.pop()                
        for Mdecay in Out_BRs:
            if Mdecay.br == 0.: continue 
            if ID in Mdecay.fstateIDs:
                newDecays = Mdecay.compress(ID,earlyDecays[ID])
                for newdecay in newDecays: newDecayList.addDecay(newdecay)            
            else:
                newDecayList.addDecay(Mdecay)
            
        Out_BRs = newDecayList
        unstableIDs = set(Out_BRs.getAllFinalStates()).intersection(set(earlyDecays.keys()))

        
    return Out_BRs

@memoize
def axionBRs(T):
    """Compute the axion branching ratios BR(a -> ids).
    Allows for temperature dependent BRs. 
    """
       
    return DecayList()

# @memoize
def axinoBRs(T):
    """Compute the axino branching ratios BR(axino -> ids).
    Allows for temperature dependent BRs. 
    """
    if modelParameters.useRGE: Tval = T
    else: Tval = susyQ
    allBRs = axinoDecays.getAxinoBRs(Tval)
    return compressBRs(allBRs,Tval)
    
@memoize
def saxionBRs(T):
    """Compute the saxion branching ratios BR(saxion -> ids).
    Allows for temperature dependent BRs.
    """

    if modelParameters.useRGE: Tval = T
    else: Tval = susyQ
    allBRs = saxionDecays.getSaxionBRs(Tval)
    return compressBRs(allBRs,Tval)
 
@memoize   
def gravitinoBRs(T):
    """Compute the gravitino branching ratios BR(gravitino -> ids).
    Allows for temperature dependent BRs.
    """
    
    if modelParameters.useRGE: Tval = T
    else: Tval = susyQ
    allBRs = gravitinoDecays.getGravitinoBRs(Tval)
    return compressBRs(allBRs,Tval)
    
@memoize   
def neutralinoBRs(T):
    """Compute the neutralino branching ratios BR(neutralino -> ids).
    Allows for temperature dependent BRs.
    """
    
    if modelParameters.useRGE: Tval = T
    else: Tval = susyQ    
    allBRs = n1Decays.getN1BRs(Tval)    
    return compressBRs(allBRs,Tval)
    

def axionCOBRs(T):    
    return axionBRs(T)
def saxionCOBRs(T):    
    return saxionBRs(T)

