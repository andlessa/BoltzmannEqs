#!/usr/bin/env python

"""

.. module:: n1Decays
    :synopsis: Computes the neutralino1 BRs and width (if unstable) 

:synopsis: Computes the lightest neutralino decays
:author: Andre Lessa <lessa.a.p@gmail.com>, Hasan Baris

"""

import modelParameters
from AuxFuncs import PSlamb, memoize,getParsAt
from AuxDecays import Decay, DecayList
from parameters import Pi,mass_Z,sw2,mtau
from math import sqrt,atan,log
from scipy import integrate
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@memoize
def getN1BRs(T):
    """Compute neutralino1 BRs at temperature T, given the parameters defined in modelParameters.
    Assumes the neutralino is the MSSM LSP (but might be heavier than axino).
    Assumes all neutralino decay energy goes into radiation."""
    
    mass_axino = modelParameters.mass_axino
    Masses = modelParameters.Masses
    slhadata = modelParameters.slhadata
    modelType = modelParameters.modelType
    caYY = modelParameters.caYY
    fa = modelParameters.fa
    
    mn1 = Masses[1000022]
    maxino = abs(mass_axino)    
    if mn1 > maxino:
        if modelType == 'DFSZ':
            logger.error('Neutralino decays are not implemented for the DFSZ model')
            return False
        else:
            logger.warning('Neutralino decays for the KSVZ model are incomplete (no Wino coulings to axino)')
    else: 
        width = 0.
        decays = DecayList()
        decays.width = width
        decays.Xfraction = 1.    
        return decays
    
    
    BRs = {}  
    n1LL = n1QQ = n1NN = n1ZA = n1GaA = 0.
    
#Get relevant parameters      
    pars = getParsAt(T)
    gSU2 = pars['gSU2']
    alphaE = gSU2**2*sw2/(4.*Pi)    
    alphaY = alphaE/(1.-sw2) 
    v41 = slhadata["NMIX"][1,1]
    gammaz = 2.4952    #Z width
    sgz = sqrt(4.*Pi*alphaE/(sw2*(1.-sw2)))
    gz = (alphaY/(16.*Pi))*caYY*v41*sqrt(sw2)/fa
    gg = (alphaY/(16.*Pi))*caYY*v41*sqrt(1.-sw2)/fa

#...SM particle dependent constants (u,d,c,s,b,e,mu,tau,neutrinos)
# Fermion masses (only used for IR cut-off)
#     mfv = [0.28,0.28,4.,1.,10.,0.5*10.**(-3),105.6*10.**(-3),mtau,1*10.**(-5)]
    mfv = [0.0002,0.0002,0.0002,0.0002,4.7,0.5*10.**(-3),105.6*10.**(-3),mtau,10.**(-5)]
#Fermion charges       
    q = [2./3.,-1./3.,2./3.,-1./3.,-1./3.,-1.,-1.,-1.,0.]
#Fermion SU2 vectorial couplings    
    gvv = [1./4.-2*sw2/3,-1./4.+sw2/3,1./4.-2*sw2/3,-1./4.+sw2/3,-1./4.+sw2/3.,-1./4.+sw2,-1./4.+sw2,-1./4.+sw2,1./4.]
#Fermion SU2 axial couplings    
    gav = [-1./4.,1./4.,-1./4.,1./4.,1./4.,1./4.,1./4.,1./4.,-1./4.]
    

#neutralino -> axino + gamma:        
    n1GaA = 2.*gg**2*((mn1**2 - maxino**2)**3)/(Pi*mn1**3)
#neutralino -> axino + Z:
    if mn1 > mass_Z + maxino:  # 2-body decay                
        PS = sqrt((-maxino+mn1-mass_Z)*(maxino+mn1-mass_Z)*(-maxino+mn1+mass_Z)*(maxino+mn1+mass_Z))
        n1ZA = gz**2*(maxino+mn1-mass_Z)*(maxino+mn1+mass_Z)*PS*(2.*((maxino-mn1)**2) + mass_Z**2)/(Pi*mn1**3)
    else:      # 3-body decays to axino-fermion-fermion (through off-shell Z plus Z/gamma interference term)
        for iF, mf in enumerate(mfv):
            qe = q[iF]*sqrt(4.*Pi*alphaE)
            gv = gvv[iF]
            ga = gav[iF]     
            x1a = 4.*mf**2/mn1**2
            x1b = 1. - maxino**2/mn1**2 - 2.*mf*maxino/mn1**2
            resZZ = integrate.romberg(n1DiffRate,x1a,x1b,args=(mn1,maxino,gammaz,mf,ga,gv,gg,gz,qe,sgz,'ZZ'),vec_func=False,rtol=0.01,divmax=30)            
            resgZ = integrate.romberg(n1DiffRate,x1a,x1b,args=(mn1,maxino,gammaz,mf,ga,gv,gg,gz,qe,sgz,'gaZ'),vec_func=False,rtol=0.01,divmax=30) 
            res = 0.004031441804*sgz**2*gz**2*mn1*resZZ 
            res += 0.002015720902*gv*sgz*gz*mn1**5*qe*gg*resgZ
            if iF <= 4: n1QQ += 3.*res      # 3 = color factor
            elif iF <= 7: n1LL += res   #leptons
            else: n1NN = 3.*res      # 3 = three neutrino flavors    
        
#...Get total width and BRs:
    Gamma = {}
    Gamma['n1-Z-A'] = n1ZA
    Gamma['n1-Ga-A'] = n1GaA
    Gamma['n1-A-lep-lep'] = n1LL
    Gamma['n1-A-nu-nu'] = n1NN
    Gamma['n1-A-q-q'] = n1QQ
    width = sum(Gamma.values())
    decays = DecayList()
    decays.width = width
    decays.Xfraction = 1.
    if width > 0.:
        decays.addDecay(Decay('neutralino',['radiation','axino'],(n1GaA+n1ZA+n1LL+n1NN+n1QQ)/width))
     
    return decays        


def n1DiffRate(x1,mn1,maxino,gammaz,mf,ga,gv,gg,gz,qe,sgz,dia=None):
    """Differential decay rate for neutralino1 -> fermion+fermion + axino, as a function of the fermion momentum"""

    mua = maxino**2/mn1**2
    muq = mf**2/mn1**2

    try:
#...logs:        
        l1 = (x1*(-1. + mua + x1) + 2.*muq*(-3. + mua + 2.*x1) -1.*sqrt(-4.*muq + x1**2)*sqrt(mua**2 + (-1. + x1)**2 
        + 2.*mua*(-1. - 2.*muq + x1)))/(-1. + muq + x1)
        
        l2 = (x1*(-1. + mua + x1) + 2.*muq*(-3. + mua + 2.*x1) + sqrt(-4.*muq + x1**2)*sqrt(mua**2 + (-1. + x1)**2 
        + 2.*mua*(-1. - 2.*muq + x1)))/(-1. + muq + x1)
        
        l3 = (2.*gammaz**2*mass_Z**2 + 2.*mass_Z**4*(-1. + x1)**2 + mn1**4*x1*(-1. + mua + x1)*(x1*(-1. + mua + x1) 
        - 1.*sqrt(-4.*muq + x1**2)*sqrt(mua**2 + (-1. + x1)**2 + 2.*mua*(-1. - 2.*muq + x1))) 
        + 2.*muq**2*(gammaz**2*mass_Z**2 + mass_Z**4 -2.*mn1**2*mass_Z**2*(-3. + mua + 2.*x1) +mn1**4*(mua**2 
        + (3. - 2.*x1)**2 + mua*(-2. + 4.*x1)) ) + 2.*mass_Z**2*(gammaz**2*(-2. + x1)*x1 +mn1**2*(-1. + x1)*(-1.*x1*(-1. 
        + mua + x1) + sqrt(-4.*muq + x1**2)*sqrt(mua**2 + (-1. + x1)**2 +2.*mua*(-1. - 2.*muq + x1)))) 
        - 2.*muq*(-2.*mass_Z**2*(gammaz**2 + mass_Z**2)*(-1. + x1) +mn1**2*mass_Z**2*(6. + mua*(-2. + 3.*x1) 
        + x1*(-11. + 5.*x1) - 1.*sqrt(-4.*muq + x1**2)*sqrt(mua**2 + (-1. + x1)**2 +2.*mua*(-1. - 2.*muq + x1))) 
        + mn1**4*((-1. + mua)**2 - 1.*(-5. + mua)*(-1. + mua)*x1 - 2.*(-3. + mua)*x1**2 - 2.*x1**3 - 3.*sqrt(-4.*muq 
        + x1**2)*sqrt(mua**2 + (-1. + x1)**2 +2.*mua*(-1. - 2.*muq + x1)) + mua*sqrt(-4.*muq + x1**2)*sqrt(mua**2 
        + (-1. + x1)**2 +2.*mua*(-1. - 2.*muq + x1)) + 2.*x1*sqrt(-4.*muq + x1**2)*sqrt(mua**2 + (-1. + x1)**2 
        + 2.*mua*(-1. - 2.*muq + x1)))))/(-1. + muq + x1)**2
        
        l4 = (2.*gammaz**2*mass_Z**2 + 2.*mass_Z**4*(-1. + x1)**2 +mn1**4*x1*(-1. + mua + x1)*(x1*(-1. + mua + x1) 
        + sqrt(-4.*muq + x1**2)*sqrt(mua**2 + (-1. + x1)**2 +2.*mua*(-1. - 2.*muq + x1))) +2.*muq**2*(gammaz**2*mass_Z**2 
        + mass_Z**4 - 2.*mn1**2*mass_Z**2*(-3. + mua + 2.*x1) + mn1**4*(mua**2 + (3. - 2.*x1)**2 +mua*(-2. + 4.*x1))) 
        + 2.*mass_Z**2*(gammaz**2*(-2. + x1)*x1 - 1.*mn1**2*(-1. + x1)*(x1*(-1. + mua + x1) +sqrt(-4.*muq 
        + x1**2)*sqrt(mua**2 + (-1. + x1)**2 + 2.*mua*(-1. - 2.*muq + x1)))) +2.*muq*(2.*mass_Z**2*(gammaz**2 
        + mass_Z**2)*(-1. + x1) - 1.*mn1**2*mass_Z**2*(6. + mua*(-2. + 3.*x1) +x1*(-11. + 5.*x1) +sqrt(-4.*muq 
        + x1**2)*sqrt(mua**2 + (-1. + x1)**2 + 2.*mua*(-1. - 2.*muq + x1))) + mn1**4*(-1. + mua**2*(-1. + x1) 
        - 3.*sqrt(-4.*muq + x1**2)*sqrt(mua**2 + (-1. + x1)**2 + 2.*mua*(-1. - 2.*muq + x1)) +mua*(2. + 2.*(-3. + x1)*x1 
        + sqrt(-4.*muq + x1**2)*sqrt(mua**2 + (-1. + x1)**2 + 2.*mua*(-1. - 2.*muq + x1))) +x1*(5. + 2.*(-3. + x1)*x1 
        + 2.*sqrt(-4.*muq + x1**2)*sqrt(mua**2 + (-1. + x1)**2 + 2.*mua*(-1. - 2.*muq + x1))))))/(-1. + muq + x1)**2
        
        l5 = l1
        l6 = l2
        
        l7 = (2.*gammaz**2*mass_Z**2 +2.*mass_Z**4*(-1. + x1)**2 +mn1**4*x1*(-1. + mua + x1)*(x1*(-1. + mua + x1) 
        - 1.*sqrt(-4.*muq + x1**2)*sqrt(mua**2 + (-1. + x1)**2 + 2.*mua*(-1. - 2.*muq + x1))) 
        + 2.*muq**2*(gammaz**2*mass_Z**2 + mass_Z**4 - 2.*mn1**2*mass_Z**2*(-3. + mua + 2.*x1) + mn1**4*(mua**2 + (3. 
        - 2.*x1)**2 + mua*(-2. + 4.*x1))) +2.*mass_Z**2*(gammaz**2*(-2. + x1)*x1 + mn1**2*(-1. + x1)*(-1.*x1*(-1. 
        + mua + x1) + sqrt(-4.*muq + x1**2)*sqrt(mua**2 + (-1. + x1)**2 + 2.*mua*(-1. - 2.*muq + x1)))) 
        - 2.*muq*(-2.*mass_Z**2*(gammaz**2 + mass_Z**2)*(-1. + x1) + mn1**2*mass_Z**2*(6. + mua*(-2. + 3.*x1) 
        + x1*(-11. + 5.*x1) - 1.*sqrt(-4.*muq + x1**2)*sqrt(mua**2 + (-1. + x1)**2 + 2.*mua*(-1. - 2.*muq + x1))) 
        + mn1**4*((-1. + mua)**2 - 1.*(-5. + mua)*(-1. + mua)*x1 - 2.*(-3. + mua)*x1**2 - 2.*x1**3 - 3.*sqrt(-4.*muq 
        + x1**2)*sqrt(mua**2 + (-1. + x1)**2 + 2.*mua*(-1. - 2.*muq + x1)) + mua*sqrt(-4.*muq + x1**2)*sqrt(mua**2 + (-1. 
        + x1)**2 + 2.*mua*(-1. - 2.*muq + x1)) + 2.*x1*sqrt(-4.*muq + x1**2)*sqrt(mua**2 + (-1. + x1)**2 
        + 2.*mua*(-1. - 2.*muq + x1)))))/ (-1. + muq + x1)**2
        
        l8 = l4
#...gamma-gamma contribution:
        dgam1 = (2.*((-1. + mua)**2 - 2.*(-1. + sqrt(mua))**2*muq +2.*(-1. + mua)*x1 
        + 2.*x1**2)*log(l1/l2) + 2.*((sqrt(-4.*muq + x1**2)*((-1. + mua)*(-2.*mua**1.5*(1. + muq) + mua**2*(1. + muq) 
        -2.*mua*(1. + muq)*(2. + muq) +(1. + muq)*(3. + 2.*muq) +2.*sqrt(mua)*(1. + 9.*muq)) + ((-1. + mua)**2*(5. 
        - 2.*sqrt(mua) + mua) +2.*(13. + 12.*sqrt(mua) - 14.*mua - 4.*mua**1.5 + mua**2)*muq)*x1 + 4.*(1. - 2.*mua**1.5 
        + mua**2 -2.*sqrt(mua)*(-1. + muq) - 7.*muq +mua*(-2. + 3.*muq))*x1**2 +2.*(-5. - 2.*sqrt(mua) + 5.*mua 
        + 4.*muq)*x1**3 + 4.*x1**4)*sqrt(mua**2 + (-1. + x1)**2 + 2.*mua*(-1. - 2.*muq + x1)))/((-1. + muq + x1)*((-1. 
        + mua)*(-1. + mua + (-9. + mua)*muq) + ((-1. + mua)**2 + 4.*(-3. + mua)*muq)*x1 + 4.*(-1. + mua + muq)*x1**2 
        + 2.*x1**3))))
#...z-z contribution:
        dgam2 = ((-2.*mn1**2*sqrt(-4.*muq + x1**2)*sqrt(mua**2 + (-1. + x1)**2 
        +2.*mua*(-1. - 2.*muq + x1))*(ga**2*(1. + 2.*sqrt(mua) - 1.*mua + 2.*muq - 2.*x1) -1.*gv**2*(-1. - 2.*sqrt(mua) 
        + mua + 2.*x1)))/ (-1. + muq + x1) + (2.*(ga**2*(2.*muq*((-1. + mua)**2*mn1**4 + (gammaz**2 
        + 4.*sqrt(mua)*mn1**2)*mass_Z**2 - 1.*mass_Z**4) +mass_Z**2*(-1.*gammaz**2*(-1. - 2.*sqrt(mua) + mua + 2.*x1) 
        + mass_Z**2*(-1. - 2.*sqrt(mua) + mua + 2.*x1) - 1.*mn1**2*((-1. + mua)**2 + 2.*(-1. + mua)*x1 + 2.*x1**2))) 
        - 1.*gv**2*(2.*(-1. + mua)**2*muq*mn1**4 -2.*(-1. + sqrt(mua))**2*muq*mn1**2*mass_Z**2 +mass_Z**2*(gammaz**2*(-1. 
        - 2.*sqrt(mua) + mua + 2.*x1) -1.*mass_Z**2*(-1. - 2.*sqrt(mua) + mua + 2.*x1) + mn1**2*((-1. + mua)**2 + 2.*(-1. 
        + mua)*x1 +2.*x1**2))))*acot((2.*gammaz*mass_Z*(-1. + muq + x1))/(2.*mass_Z**2*(-1. + x1) - 1.*mn1**2*x1*(-1. 
        + mua + x1) + mn1**2*sqrt(-4.*muq + x1**2)*sqrt(mua**2 + (-1. + x1)**2 +2.*mua*(-1. - 2.*muq + x1)) 
        +2.*muq*(mass_Z**2 - 1.*mn1**2*(-3. + mua + 2.*x1)))) +2.*(gv**2*(2.*(-1. + mua)**2*muq*mn1**4 -2.*(-1. 
        + sqrt(mua))**2*muq*mn1**2*mass_Z**2 +mass_Z**2*(gammaz**2*(-1. - 2.*sqrt(mua) + mua + 2.*x1) - 1.*mass_Z**2*(-1. 
        - 2.*sqrt(mua) + mua + 2.*x1) + mn1**2*((-1. + mua)**2 + 2.*(-1. + mua)*x1 +2.*x1**2))) + ga**2*(2.*muq*(-1.*(-1. 
        + mua)**2*mn1**4 - 1.*(gammaz**2 + 4.*sqrt(mua)*mn1**2)*mass_Z**2 + mass_Z**4) +mass_Z**2*(gammaz**2*(-1. 
        - 2.*sqrt(mua) + mua + 2.*x1) - 1.*mass_Z**2*(-1. - 2.*sqrt(mua) + mua + 2.*x1) + mn1**2*((-1. + mua)**2 
        + 2.*(-1. + mua)*x1 +2.*x1**2))))*acot((2.*gammaz*mass_Z*(-1. + muq + x1))/(2.*mass_Z**2*(-1. + x1) 
        + 2.*muq*(mass_Z**2 - 1.*mn1**2*(-3. + mua + 2.*x1)) - 1.*mn1**2*(x1*(-1. + mua + x1) +sqrt(-4.*muq 
        + x1**2)*sqrt(mua**2 + (-1. + x1)**2 + 2.*mua*(-1. - 2.*muq + x1))))) +gammaz*mass_Z*(ga**2*(-2.*mass_Z**2*(-1. 
        - 2.*sqrt(mua) + mua - 2.*muq + 2.*x1) +mn1**2*(1. + mua**2 - 8.*sqrt(mua)*muq + 2.*mua*(-1. + x1) + 2.*(-1. 
        + x1)*x1)) +gv**2*(-2.*mass_Z**2*(-1. - 2.*sqrt(mua) + mua + 2.*x1) +mn1**2*((-1. + mua)**2 - 2.*(-1. 
        + sqrt(mua))**2*muq + 2.*(-1. + mua)*x1 + 2.*x1**2)))*(log(l3/l4)))/(gammaz*mass_Z))
#...z-gamma contribution:
        dgam12 = ((4.*(-1. - 2.*sqrt(mua) + mua + 2.*x1)*sqrt(-4.*muq 
        + x1**2)*sqrt(mua**2 + (-1. + x1)**2 + 2.*mua*(-1. - 2.*muq + x1)) )/(mn1**2*(-1. + muq + x1)) 
        + (2.*(2.*gammaz*(2.*(-1. + mua)**2*muq*mn1**4 + mass_Z**2*(gammaz**2 + mass_Z**2)*(-1. - 2.*sqrt(mua) + mua 
        + 2.*x1))*atan((2.*gammaz*mass_Z*(-1. + muq + x1))/(2.*mass_Z**2*(-1. + x1) - 1.*mn1**2*x1*(-1. + mua + x1) 
        + mn1**2*sqrt(-4.*muq + x1**2)*sqrt(mua**2 + (-1. + x1)**2 +2.*mua*(-1. - 2.*muq + x1)) +2.*muq*(mass_Z**2 
        - 1.*mn1**2*(-3. + mua + 2.*x1)))) +2.*gammaz*(-2.*(-1. + mua)**2*muq*mn1**4 - 1.*mass_Z**2*(gammaz**2 
        + mass_Z**2)*(-1. - 2.*sqrt(mua) + mua + 2.*x1))*atan((2.*gammaz*mass_Z*(-1. + muq + x1))/(2.*mass_Z**2*(-1. 
        + x1) + 2.*muq*(mass_Z**2 - 1.*mn1**2*(-3. + mua + 2.*x1)) - 1.*mn1**2*(x1*(-1. + mua + x1) +sqrt(-4.*muq 
        + x1**2)*sqrt(mua**2 + (-1. + x1)**2 + 2.*mua*(-1. - 2.*muq + x1))))) +mass_Z*(-4.*(-1. 
        + mua)**2*muq*mn1**4*log(l5/l6) + (2.*(-1. + mua)**2*muq*mn1**4 -2.*(-1. + sqrt(mua))**2*muq*mn1**2*mass_Z**2 
        + mass_Z**2*(-1.*mass_Z**2*(-1. - 2.*sqrt(mua) + mua + 2.*x1) + mn1**2*((-1. + mua)**2 + 2.*(-1. + mua)*x1 
        + 2.*x1**2)) +gammaz**2*(-1.*mass_Z**2*(-1. - 2.*sqrt(mua) + mua + 2.*x1) + mn1**2*((-1. + mua)**2 -2.*(-1. 
        + sqrt(mua))**2*muq +2.*(-1. + mua)*x1 + 2.*x1**2)))*(log(l7/l8)))))/ (mn1**4*mass_Z*(gammaz**2 + mass_Z**2)))
        
        if dia == 'gaga': return dgam1
        elif dia == 'ZZ': return dgam2
        elif dia == 'gaZ': return 2.*dgam12
#         return dgam1 + dgam2 + 2.*dgam12
    except: return 0.

def acot(x):
    
    return atan(1./x)
