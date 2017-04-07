#!/usr/bin/env python

"""

.. module:: modelParameters
    :synopsis: This module defines all relevant input parameters for the model and some methods to obtain \
    other relevant information (SUSY spectrum and parameters, neutralino sigma.v, etc)
    
:synopsis: This module defines all relevant input parameters for the model and some methods to obtain \
    other relevant information (SUSY spectrum and parameters, neutralino sigma.v, etc)    
:author: Andre Lessa <lessa.a.p@gmail.com>

"""

import subprocess, tempfile
from scipy import interpolate
import pyslha,os
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


MSSMparticles={
1000021 : "gluino", 1000022: "neutralino", 1000023 : "neutralino", 1000025 : "neutralino", 1000035 : "neutralino", 
1000024 : "chargino", 1000037 : "chargino", 1000039 : "gravitino", 1000001 : "squark", 1000002 : "squark", 
1000003 : "squark", 1000004 : "squark", 2000001 : "squark", 2000002 : "squark", 2000003 : "squark", 
2000004 : "squark", 1000005 : "squark", 2000005 : "squark", 1000006 : "squark", 2000006 : "squark", 
1000011 : "slepton", 1000013 : "slepton", 1000015 : "slepton", 2000011 : "slepton", 2000013 : "slepton", 
2000015 : "slepton", 1000012 : "sneutrino", 1000014 : "sneutrino", 1000016 : "sneutrino", 
2000012 : "sneutrino", 2000014 : "sneutrino", 2000016 : "sneutrino", -1000021 : "gluino",
25 : "higgs0", 35 : "higgs0", 36 : "higgs0", 37 : "higgs+" 
}


useRGE = True #If True, use the RGE values from rgedata to compute parameters at T. If False, use values at the SUSY scale
   
def getIsajetData(codeDir,sigmav_datafile,slhafile,rge_datafile,isajet_inputfile=None):
    """
    Obtains additional model specific data.
    
    :param codeDir: location of the code folder (where getIsajetData.x is)
    :param sigmav_datafile: File containing T,sigma.v(T) points to be interpolated and compute the neutralino annihilation cross-section
    :param slhafile: File containing the MSSM spectrum and parameters
    :param rge_datafile: File containing T,M1(T),M2(T),M3(T),g'(T),gSU2(T),gSU3(T)
    :param isajet_inputfile: File containing isasugra type input to run isajet, isared and compute the relevant data
    
    If sigmav_datafile, slhafile and/or rge_datafile does not exist, run isajet and isared with isajet_input
    and compute the missing data. Save the missing data to the respective file.
    If sigmav_datafile, slhafile and/or rge_datafile = None, run isajet and isared with isajet_input
    but does not save the missing data to file.    
    """
    
    global Masses,slhadata,sigmaFunc,rgeFunc
    
    tempfiles = []
    datafiles = {'sigmavfile' : sigmav_datafile,'slhafile' : slhafile, 'rgefile' : rge_datafile}
    runIsajet = False
    for f in datafiles.values():
        if not f or not os.path.isfile(f):
            if not isajet_inputfile or not os.path.isfile(isajet_inputfile):
                logger.error("One of the data files is missing and the isajet input file was not defined.")
                return False
            runIsajet = True
            
#Run isajet/isared if needed:
    if runIsajet:
        if not codeDir or not os.path.isdir(codeDir) or not os.path.isfile(codeDir+'/getIsajetData.x'):
            logger.error("Could not find getIsajetData.x in "+codeDir)
            return False
        isa_file = open(isajet_inputfile, 'r')
        isa_in = isa_file.read()
        isa_file.close()
        isa_head = ""        
        for tag,f in datafiles.items():
            if not f:  #Do not save data to file (store in a temp file)
                fname = gettempfile()
                tempfiles.append(fname)
                datafiles[tag] = fname                
            elif not os.path.isfile(f): fname = f  #Save data to file rgedata
            else: fname = None   #Use existing datafile
            isa_head += tag+","+str(fname)+"\n"    

#Create getIsajetData inputfile with the names of the files to be used (if any) and the isasugra input:
        isa_in = isa_head + isa_in
        tmp = gettempfile()
        tempfiles.append(tmp)
        isa_file = open(tmp,'w')
        isa_file.write(isa_in)
        isa_file.close()
        data_in = gettempfile()
        tempfiles.append(data_in)
        f = open(data_in,'w')
        f.write("\'"+tmp+"\'"+'\n')
        f.close()
#Run getIsajetData.x:
        subprocess.call(codeDir+"getIsajetData.x < "+data_in,shell=True)

#Store SLHA data:
    data = pyslha.readSLHAFile(datafiles['slhafile'])      
    Masses = {}
    for pid,mass in data[0]['MASS'].items():        
        if pid and mass: Masses[pid] = abs(mass)
    slhadata = data[0]

#Get sigmaV data:
    x = "T"
    y = ["sigmaV"]
    sigmaVdata = getDataFrom(datafiles['sigmavfile'],x,y)
#Create sigmaV function
    sigF = interpolate.interp1d(sigmaVdata[x],sigmaVdata[y[0]],kind='linear') 
    def sigmaFunc(T):
        if T > sigmaVdata["T"][-1]: return sigmaVdata["sigmaV"][-1]
        elif T < sigmaVdata["T"][0]: return sigmaVdata["sigmaV"][0]
        else: return sigF(T)
        
#Get rge data
    x = "Q"
    y = ["M1","M2","M3","gPr","gSU2","gSU3"]
    rgedata = getDataFrom(datafiles['rgefile'],x,y)        
#Create rge function
    rgeF = interpolate.interp1d(rgedata[x],[rgedata[yl] for yl in y],kind='linear')
    def rgeFunc(T):
        res = {}
        if T > rgedata["Q"][-1]:
            for yl in y: res[yl] = rgedata[yl][-1]
        elif T < rgedata["Q"][0]:
            for yl in y: res[yl] = rgedata[yl][0]
        else:
            r = rgeF(T)            
            for i,yl in enumerate(y): res[yl] = r[i]
        return res
        
#Cleanup temp files:
    for temp in tempfiles: os.remove(temp)
        
    return True


def gettempfile():
    """Simple method to create a tempfile in the current folder. Returns the tempfile name."""
    
    
    tmp = tempfile.mkstemp(dir=os.getcwd())
    os.close(tmp[0])
    return tmp[1]

def getDataFrom(fname,x,y,sort=True):
    """Reads fname and returns a dictionary with a list of values for each entry in y. The
    number of labels in x+y  must match the number of columns in the file"""
    
           
    f = open(fname,'r')
    data = f.read()
    f.close()
    data  = data.split('\n')
    datapts = []    
    for data_pt in data:
        error = False
        if not data_pt: continue
        pt = data_pt.split()        
        pt = [eval(p) for p in pt]                
        if len(pt) != len(y)+1: error = True
        for p in pt:            
            if type(p) != type(float()): error = True
        if error:
            logger.error("Error reading "+fname+"  data.")
            return False
        else:            
            datapts.append(pt)
            
#Sort sigmaVdata by temperature
    datapts = sorted(datapts)
    dataDic = {}    
    dataDic[x] = []
    for ylabel in y: dataDic[ylabel] = []
    for pt in datapts:
        dataDic[x].append(pt[0])    
        for i,ylabel in enumerate(y): dataDic[ylabel].append(pt[i+1])
        
    return dataDic
