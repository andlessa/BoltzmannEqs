#!/usr/bin/env python

"""

.. module:: AuxFuncs
    :synopsis: This module provides auxiliary functions 

:synopsis: This module provides auxiliary functions
:author: Andre Lessa <lessa.a.p@gmail.com>

"""

from scipy import optimize, integrate, interpolate
from math import exp, sqrt, log10, log
from parameters import Pi, MP
import modelParameters
import logging
from functools import wraps
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
import warnings

warnings.filterwarnings('error')

def memoize(func):
    """ A wrapping to store in cache the results of a function, for expensive functions.
    Ignores all but the first argument of the function."""
    
    cache = {}
    @wraps(func)
    def wrap(*args):
        if not str(args[0]) in cache: cache[str(args[0])] = func(*args)
        return cache[str(args[0])]
    return wrap


def plainmemoize(func):
    """ A wrapping to store in cache the results of a function, for expensive functions.
    Ignores all the inputs of the function and returns always the same result."""
    
    cache = {}
    @wraps(func)
    def wrap(*args):
        if not cache: cache[str(args[0])] = func(*args)
        return cache.values()[0]
    return wrap

def printSummary(compList,TF):
    """Prints basic summary of solutions."""
#Solution summary:
    print '----Summary----'
    print 'TF=',TF
    for comp in compList:         
        mindelta, Tlast = None, TF
        if comp.Tdecay: Tlast = max(comp.Tdecay,TF)    
        for iT,T in enumerate(comp.evolveVars['T']):   #Get point closest to Tlast        
            if abs(T-Tlast) < mindelta or mindelta is None: iF, mindelta = iT, abs(T-Tlast)
        rhoF = comp.evolveVars['rho'][iF]
        nF = comp.evolveVars['n'][iF]
        Tfinal = comp.evolveVars['T'][iF]        
        omega = getOmega(comp,rhoF,nF,Tfinal)        
        if not comp.Tdecay: tag = '(@TF)'
        else: tag = '(@decay)'
        print comp.label+':','T(osc)=',comp.Tosc,' | T(decay)=',comp.Tdecay,' | Omega h^2',tag,'=',omega
        print '-'*len(comp.label)
    
    print 'Delta Neff (@TF) =',sum([getDNeff(comp,TF) for comp in compList])



def getParsAt(T):
    """Try to obtain soft parameters and gauge couplings at temperature T from
    the RGE function defined in modelParameters.rgeFunc.
    If the RGE function does not contain the running of a parameter, use its (default) SUSY scale
    value from the SLHA input data"""
    
    slhadata = modelParameters.slhadata
    allPars = {'gPr' : slhadata["GAUGE"][1],'gSU2' : slhadata["GAUGE"][2],'gSU3' : slhadata["GAUGE"][3]
               ,'M1' : abs(slhadata["MSOFT"][1]),'M2' : abs(slhadata["MSOFT"][2])
               ,'M3' : abs(slhadata["MSOFT"][3]),'mu' : slhadata["HMIX"][1], 'tanb' :  slhadata["HMIX"][2]}

    if not modelParameters.useRGE: return allPars    
    rgePars = modelParameters.rgeFunc(T)
    parsT = {}
    for par in allPars:
        if par in rgePars: parsT[par] = rgePars[par]
        else: parsT[par] = allPars[par]
        
    return parsT



def getOmega(comp,rho,n,T):
    """Compute relic density today, given the component, number density and energy density\
    at temperature T. """
    
    if comp.Tdecay and comp.Tdecay > T: return 0.
    
    Ttoday = 2.3697*10**(-13)*2.725/2.75  #Temperature today
    rhoh2 = 8.0992*10.**(-47)   # value of rho critic divided by h^2
    dx = (1./3.)*log(gSTARS(T)/gSTARS(Ttoday)) + log(T/Ttoday)   #dx = log(R/R_today), where R is the scale factor
    nToday = n*exp(-3.*dx)
    ns = log((2*Pi**2/45)*T**3)  #entropy (constant) 
    
    if comp.Type == 'CO': return nToday*comp.mass(Ttoday)/rhoh2  #CO components have trivial (non-relativistic) solution 
               
    R0 = rho/n    
    Rmin = R0*exp(-dx)    #Minimum value for rho/n(Ttoday) (happens if component is relativistic today)
    Pmin = getPressure(comp.mass(Ttoday),Rmin*nToday,nToday)
      
           
    if abs(Pmin - Rmin*nToday/3.)/(Rmin*nToday/3.) < 0.01: RToday = Rmin  #Relativistic component today
    else:
        def Rfunc(R,x):            
            TF = getTemperature(x,ns)
            nF = n*exp(-3*x)   #Number density at x (decoupled solution)
            rhoF = R*nF         #Energy density at x                        
            return -3*getPressure(comp.mass(TF),rhoF,nF)/nF
        RToday = integrate.odeint(Rfunc, R0, [0.,24.], atol = comp.mass(Ttoday)/10.)[1][0]  #Solve decoupled ODE for R=rho/n
   
    return RToday*nToday/rhoh2


def getDNeff(comp,TF):
    """Computes the contribution from component comp to the number of effective neutrinos at temperature TF.
    Can only be used after the Boltzmann equations have been solved and the solutions stored in comp.evolveVars.
    Gives zero if T > 1 MeV (where neutrinos are still coupled)."""

#Ignore component if it has decayed before TF        
    if comp.Tdecay and comp.Tdecay > TF: return 0.
#Get the number and energy densities of comp at T:
    mindelta = TF    
    for iT,T in enumerate(comp.evolveVars['T']):
        if abs(T-TF) < mindelta: iF,mindelta = iT,abs(T-TF)
    rho = comp.evolveVars['rho'][iF]
    n = comp.evolveVars['n'][iF]
    T = comp.evolveVars['T'][iF]
    mass = comp.mass(T)
    if T > 10.**(-3): return 0.    
    if mass == 0. or (n and rho and rho/(n*mass) > 2.): rhoRel = rho    
    else: rhoRel = 0.
    DNeff = rhoRel/(((Pi**2)/15)*(7./8.)*((4./11.)**(4./3.))*T**4)
    return DNeff


def Hfunc(x, nv, rhov, NS, sw):
    """Compute the Hubble parameter, given the variables x=log(R/R0) and ni, rhoi and NS=log(S/S0) """
    
    T = getTemperature(x, NS)
    
    rhoActive = []
    for i, rho in enumerate(rhov):
        if not sw[i]: continue
        rhoActive.append(rho)  # energy density of each active component     
    rhoRad = (Pi**2/30)*gSTAR(T)*T** 4  # thermal bath's energy density    
    rhoActive.append(rhoRad)
    rhoTot = sum(rhoActive)  # Total energy density    
    H = sqrt(8*Pi*rhoTot/3)/MP
    
    return H


@memoize
def getTemperature(x,NS):
    """Computes the temperature for the thermal bath from xeff = NS - 3*x, where
    x = log(R) and NS = log(S)."""

    xeff = NS - 3*x
    
    def Tfunc(T):
        """Auxiliary function, its zero implicitly determines T"""        
        return ((2.*Pi**2)/45)*gSTARS(T)*T**3 - exp(xeff)

    Tmax = (45./(2.*Pi**2))**(1./3.)*exp(xeff/3)  # Upper limit on T (since gSTARS > 1 always)
    Tmin = (45./(2.*Pi**2*230.))**(1./3.)*exp(xeff/3)  # Lower limit on T (since gSTARS < 225 always)
    
    if gSTARS(Tmin) == gSTARS(Tmax): return (((2.*Pi**2)/45)*gSTARS(Tmin))**(-1./3.)*exp(xeff/3)
        
    return optimize.brenth(Tfunc, Tmin, Tmax)


def getPressure(mass, rho, n):
    """Computes the pressure for a component, given its mass, its energy density and its number density"""

    R = rho/n    
    if R > 11.5*mass: return n*(R/3)  # Ultra relativistic limit
    if R <= mass: return 0.  # Ultra non-relativistic limit
    
# Expansion coefficients for relativistic/non-relativistic transition    
    aV = [-0.345998, 0.234319, -0.0953434, 0.023657, -0.00360707, 0.000329645, -0.0000165549, 3.51085*10.**(-7)]    
    Prel = n*(R/3)  # Relativistic pressure
    Pnonrel = (2.*mass/3.)*(R/mass - 1.)  # Non-relativistic pressure
    Pnonrel += mass*sum([ai*(R/mass - 1.)**(i+2)  for i, ai in enumerate(aV)])
    Pnonrel *= n
        
    return min(Prel, Pnonrel)  # If P2 > P1, it means ultra relativistic limit applies -> use P1
    
def gSTAR(T):
    """Computes the number of relativistic degrees of freedom in thermal equilibrium with the thermal\
    bath, including the full MSSM spectrum, except for the lightest neutralino.\
    To save time, for the first call generates exact data points describing well the function
    and store them in gSTARdata. For all subsequent calls, interpolate between the data points.
    """
    
    global gSTARdata, interpF

# If data already exists, use it    
    try:        
        gSTARdata          
# If not, generate data points    
    except NameError:
        try:
            Masses = modelParameters.Masses
        except:
            Masses = None
        if not Masses:
            logger.error("A MSSM spectrum must be defined in modelParameters")
            return False
        
        gSTARdata = []
        maxMass = max(Masses.values())  # Use maximum mass value to get maximum relevant temperature
        Tmax = 10.*maxMass  # Get maximum T value (gSTAR is constant above this values)
        Tmin = 10.**(-5)  # gSTAR is constant below 10 keV
        for i in range(31):
            T = 10.**(log10(Tmin) + log10(Tmax / Tmin) * (i / 30.))
            gSTARdata.append([T, gSTARexact(T)])
        done = False
# Generate enough points until the linear interpolation of points describes gSTAR with a better than 10% accuracy     
        while not done:
            done = True
            Tpts = [pt[0] for pt in gSTARdata]
            gSTARpts = [pt[1] for pt in gSTARdata]
            gFunc = interpolate.interp1d(Tpts, gSTARpts, kind='linear')
            for ipt, pt in enumerate(Tpts[:-1]):                
                Tnew = (Tpts[ipt] + Tpts[ipt + 1]) / 2
                gnew = gSTARexact(Tnew)
                if abs(gnew - gFunc(Tnew)) / gnew > 0.1:  # If the interpolation is not good, add the point and try again
                    gSTARdata.insert(ipt + 1, [Tnew, gnew])
                    done = False
            gSTARdata = sorted(gSTARdata)
            Tpts = [pt[0] for pt in gSTARdata]
            gSTARpts = [pt[1] for pt in gSTARdata]
            interpF = interpolate.interp1d(Tpts, gSTARpts, kind='linear')

    if T >= gSTARdata[-1][0]: return gSTARdata[-1][1]  # gSTAR should be constant above Tmax and below Tmin
    elif T <= gSTARdata[0][0]: return gSTARdata[0][1]
    else: return float(interpF(T))



@memoize
def gSTARexact(T, interpol=True):
    """Computes exactly the number of relativistic degrees of freedom in thermal equilibrium with the thermal\
    bath, including the full MSSM spectrum, except for the lightest neutralino.\
    interpol turns on/off the interpolation around the QCD phase trasition region."""

    def gstarFunc(x, dof):
        """Auxiliary function to compute the contribution from a single particle to gSTAR.\
        x = mass/T, dof = number of degrees of freedom (positive/negative for bosons/fermions)"""
                
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

    Masses = modelParameters.Masses
    MSSMparticles = modelParameters.MSSMparticles
    
# Define fermi/bose statistics for each state and number of degrees of freedom:
    DoF = {"squark" : 6, "slepton" : 2, "sneutrino" : 1, "gluino" :-16, "neutralino" :-2, "chargino" :-4,
           "higgs0" : 1, "higgs+" : 2}
    gstar = 0.
# Add up MSSM degrees of freedom (except for neutralino1)    
    for pid in Masses:        
        if pid == 1000022 or pid == 1000039: continue # Skip neutralino1 and gravitino 
        if not pid in MSSMparticles: continue  # Skip SM particles (except Higgs)
        spart = MSSMparticles[pid]  # Sparticle name
        gstar += gstarFunc(Masses[pid] / T, DoF[spart])

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
    MassesSM = dict(MassesGauge.items() + MassesLeptons.items() + MassesHadrons.items())
    DoFSM = dict(DoFGauge.items() + DoFLeptons.items() + DoFHadrons.items())
# Add up SM degrees of freedom     
    for part in MassesSM: gstar += gstarFunc(MassesSM[part] / T, DoFSM[part])

# Correct for neutrino decoupling:
    if T <= 5.*10.**(-4):
        gstar += (-1. + (4. / 11.) ** (4. / 3.)) * gstarFunc(MassesSM["neutrino"] / T, DoFSM["neutrino"]) 
     
# Smooth discontinuous transitions:
    if interpol:
# QCD phase transition:
        finter = None
        if T < 0.3 and T > 0.15 and interpol:
            Tpts = [0.15, 0.3]
            gpts = [gSTARexact(Tpt, False) for Tpt in Tpts]        
            finter = interpolate.interp1d(Tpts, gpts, kind='linear')            
# Neutrino decoupling            
        elif T < 6.*10.**(-4) and T > 2.*10.**(-4):
            Tpts = [2.*10.**(-4), 6.*10 ** (-4)]
            gpts = [gSTARexact(Tpt, False) for Tpt in Tpts]        
            finter = interpolate.interp1d(Tpts, gpts, kind='linear')
# Replace gstar value by interpolation:
        if finter: gstar = finter(T)
    
    return gstar
    

def gSTARS(T):
    """Computes the number of relativistic degrees of freedom for computing the entropy density,\
    including the full MSSM spectrum, except for the lightest neutralino."""
    
    if T >= 10.**(-3): return gSTAR(T)
    else:  # Correct for neutrino decoupling:        
        return gSTAR(T) - (7. / 8.) * 6.*(4. / 11.) ** (4. / 3.) + (7. / 8.) * 6.*(4. / 11.)
    
def PSlamb(x,y,z):
    """Simple 3-body phase factor auxiliary function """
    
    return x**2 + y**2 + z**2 - 2.*x*y - 2.*x*z - 2.*z*y