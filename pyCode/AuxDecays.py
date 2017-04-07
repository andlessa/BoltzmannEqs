"""

.. module:: AuxDecays
    :synopsis: This module provides auxiliary classes and methods for dealing with decay data 

:synopsis: This module provides decay classes and methods
:author: Andre Lessa <lessa.a.p@gmail.com>

"""

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
import warnings
import copy

warnings.filterwarnings('error')


class Decay(object):
    """Main class to store information about a specific decay.
    
    :param motherID: ID for the decaying mother
    :param fstateIDs: list of IDs for the daughters
    :param br: branching ratio value
    
    """
    
    def __init__(self,instate=None,fstates=None,br=0.):
        self.motherID = instate
        self.fstateIDs = fstates
        self.br = br
        if self.fstateIDs: self.fstateIDs.sort()
        
    def __eq__(self,other):
                
        if type(self) != type(other): return False
        if self.motherID != other.motherID: return False
        if self.fstateIDs != other.fstateIDs: return False
        return True
    
    def __str__(self):
        return "BR("+self.motherID + " -> " + str(self.fstateIDs)+") = "+str(self.br)
    
    def compress(self,ID,IDecays):
        """Generates a new DecayList() from compressing the ID decays"""
        
        newList = DecayList()
        for dec in IDecays:
            newDecay = copy.deepcopy(self)
            newDecay.fstateIDs.remove(ID)
            newDecay.fstateIDs += dec.fstateIDs
            newDecay.br *= dec.br
            newList.addDecay(newDecay)
        
        return newList
        
class DecayList(object):
    """Main class to store information about all the decays of a specific species.
    
    :param dlist: List of Decay objects
    :param width: Total width
    :param Xfraction: fraction of energy injected in the thermal bath
    
    """
    
    def __init__(self):
        self.dlist = []
        self.width = 0.
        self.Xfraction = 0.
       
    
    def __iter__(self):
        return iter(self.dlist)
        
    def __getitem__(self, index):
        return self.dlist[index]
        
    def __setitem__(self, index, decay):
        if type(decay) != type(Decay()):
            logger.error("Input object must be a Decay() object")
            return False
        else:
            self.dlist[index] = decay
        
    def __len__(self):
        return len(self.dlist)
    
    def __eq__(self,other):    
        return self.dlist == other.dlist
    
    def __add__(self,other):
        for dec in other.dlist: self.addDecay(dec)
        
    def __str__(self):
        
        strList = "Width = "+str(self.width)+" X fraction = "+str(self.Xfraction)+"\n"
        for dec in self.dlist:
            strList += str(dec)+"\n"
        return strList
    
    def addDecay(self,decay):
       
        if type(decay) != type(Decay()): return False               
        for dec in self.dlist:
            if dec == decay:
                dec.br *= decay.br
                return        
        self.dlist.append(decay)
        return True
    
    def getAllFinalStates(self):
        
        allstates = []
        for dec in self.dlist: allstates += dec.fstateIDs
        return set(allstates)
