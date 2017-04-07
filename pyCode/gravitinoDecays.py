#!/usr/bin/env python

"""

.. module:: gravitinoDecays
    :synopsis: Computes gravitino BRs and width 

:synopsis: Computes the gravitino decays
:author: Andre Lessa <lessa.a.p@gmail.com>, Hasan Baris

"""

import modelParameters
from AuxFuncs import memoize, getParsAt
from AuxDecays import Decay, DecayList
from parameters import Pi,mass_Z,mass_W,vEW,mt,mb,mtau,MPL
from math import atan,sqrt,sin,cos
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@memoize
def getGravitinoBRs(T):
    """Compute gravitino BRs at temperature T, given the parameters defined in modelParameters.\
    Assumes all gravitino decay energy goes into radiation."""

    mass_gravitino = modelParameters.mass_gravitino
    mass_axino = modelParameters.mass_axino
    mass_saxion = modelParameters.mass_saxion
    Masses = modelParameters.Masses
    slhadata = modelParameters.slhadata

    BRs = {}
    
#Get relevant parameters  
    mgt = abs(mass_gravitino)    
    pars = getParsAt(T)
    gSU2 = pars['gSU2']
    gU1 = pars['gPr']/sqrt(3./5.)
    gz = sqrt(gU1**2 + gSU2**2)
    tanb = pars['tanb']
    beta = atan(tanb)   
    mh = Masses[25]
    mH = Masses[35]
    ma = Masses[36]
    mHp = Masses[37]

    pidsN = [1000022,1000023,1000025,1000035]
    pidsC = [1000024,1000037]
    u0st = []
    for iN in range(4):
        u0st.append([])
        for jN in range(4):
            uij = slhadata["NMIX"][iN+1,jN+1]
            if jN+1 >= 3: uij *= -1       
            u0st[iN].append(uij)
    umst,upls = [],[]        
    for iC in range(2):
        umst.append([])
        upls.append([])
        for jC in range(2):
            uij = slhadata["UMIX"][iC+1,jC+1]*1j
            vij = slhadata["VMIX"][iC+1,jC+1]*1j
            if iC+1 == 1: uij *= -1
            if jC+1 == 2: vij *= -1            
            umst[iC].append(uij)
            upls[iC].append(vij)
    utp = [[slhadata["STOPMIX"][jF+1,iF+1] for jF in range(2)] for iF in range(2)]
    ubt = [[slhadata["SBOTMIX"][jF+1,iF+1] for jF in range(2)] for iF in range(2)]
    uta = [[slhadata["STAUMIX"][jF+1,iF+1] for jF in range(2)] for iF in range(2)]
    alpha = slhadata["ALPHA"].values()[0]
    v1 = vEW*cos(beta)
    v2 = vEW*sin(beta)


    cgl, cgr, chL, chR = [[[0.]*6 for i in range(3)] for j in range(4)]
    ccl, ccr = [[[0.]*6 for i in range(4)] for j in range(2)]
    for jN in range(6):
        if jN <= 3:        
            cgl[0][jN] = cgr[0][jN] = (gSU2*u0st[jN][0]+gU1*u0st[jN][1])/gz            
            cgl[1][jN] = cgr[1][jN] = (gSU2*u0st[jN][1]-gU1*u0st[jN][0])/gz
            chL[1][jN] = gz*(-v1*u0st[jN][2] + v2*u0st[jN][3])/sqrt(2.)
            chR[1][jN] = - chL[1][jN]
            ccl[0][jN] = ccr[0][jN] = -sin(alpha)*u0st[jN][2] + cos(alpha)*u0st[jN][3]
            ccl[1][jN] = ccr[0][jN] = cos(alpha)*u0st[jN][2] + sin(alpha)*u0st[jN][3]
            ccl[2][jN] = sin(beta)*u0st[jN][2] + cos(beta)*u0st[jN][3]
            ccr[2][jN] = -ccl[2][jN]
        else:
            cgl[2][jN] =  umst[jN-4][0]
            cgr[2][jN] =  upls[jN-4][0]
            chL[2][jN] = -gSU2*v1*umst[jN-4][1]
            chR[2][jN] = gSU2*v2*upls[jN-4][1]
            ccl[3][jN] = sqrt(2.)*sin(beta)*umst[jN-4][1]
            ccr[3][jN] = sqrt(2.)*cos(beta)*upls[jN-4][1]

 
    uferms = [utp,ubt,uta]
    cfl, cfr = [[[0.]*2 for i in range(3)] for j in range(2)]
    for iF in range(2):
        for j,uferm in enumerate(uferms):
            cfl[j][iF] = sqrt(2.)*uferm[iF][0]
            cfr[j][iF] = sqrt(2.)*uferm[iF][1]
 
#...check for negative isajet masses:
    negm = []
    for pid in pidsN + pidsC:
        if slhadata["MASS"][pid] < 0: negm.append(-1.)
        else: negm.append(1.)
#...redefine left-handed couplings [ (1-gamma5)-> gamma5*(1-gamma5)=-(1-gamma5) => c*l-> -c*l ]:
    for j in range(6):
        for i in range(3):
            cgl[i][j] *= negm[j]
            chL[i][j] *= negm[j]
        for i in range(4):
            ccl[i][j] *= negm[j]
 
#...compute gravitino decay width
    Gamma = {}
    scalars = {'squarks' : [1000000+i for i in range(1,5)] + [2000000+i for i in range(1,5)], 'sbottoms' : [1000005,2000005],
            'stops' : [1000006,2000006],'sleptons' : [1000011,1000013,2000011,2000013],'staus' : [1000015,2000015],
            'sneutrinos' : [1000012,1000014,1000016]}
    fermions = {'gluino' : [1000021], 'neutralinos' : [1000022,1000023,1000025,1000035], 'charginos' : [1000024,1000037]}
    labs = dict(scalars.items() + fermions.items())
    gmult = {'gluino' : 8., 'squarks' : 3., 'sbottoms' : 3., 'stops' : 3., 'sleptons' : 1.,'staus' : 1., 'sneutrinos' : 1., 'neutralinos' : 1., 'charginos' : 1.}
    mFs = {'sbottoms' : mb, 'stops' : mt, 'staus' : mtau}
    mSs = {'neutralinos' : [mh,mH,ma,mHp], 'charginos' : [mh,mH,ma,mHp], 'gluinos' : []}
    mVs = {'neutralinos' : [0.,mass_Z,mass_W], 'charginos' : [0.,mass_Z,mass_W], 'gluinos' : [0.]}


#gravitino -> scalar + fermion decays    
    for lab,pids in labs.items():
        Gamma[lab] = 0.
        for pid in pids:
            if lab in scalars:   #Sfermion + fermion
                mscalar = Masses[pid]
                if lab in mFs: mfermion = mFs[lab]
                else: mfermion = 0.
                iF = pid/1000000 - 1
                if lab == 'stops': itype = 0
                elif lab == 'sbottoms': itype = 1
                elif lab == 'staus': itype = 2
                else: itype = None
                if not itype is None:
                    cc1 = cfl[itype][iF]**2 + cfr[itype][iF]**2
                    cc2 = 2*cfl[itype][iF]*cfr[itype][iF]                
                else:
                    cc1 = 2.
                    cc2 = 0
                Gamma[lab] += gmult[lab]*decayToFS(mgt,mfermion,mscalar,cc1,cc2)                
            elif lab in fermions:   #Higgsino + Higgs bosons
                if lab == 'gluino': continue
                iF = (fermions['neutralinos'] + fermions['charginos']).index(pid) 
                mfermion = Masses[pid]
                for j,mscalar in enumerate(mSs[lab]):
                    cc1 = ccl[j][iF]*conjg(ccl[j][iF]) + ccr[j][iF]*conjg(ccr[j][iF])
                    cc2 = ccl[j][iF]*conjg(ccr[j][iF]) + ccr[j][iF]*conjg(ccl[j][iF])
                    if mscalar == mHp: Gamma[lab] += 2*gmult[lab]*decayToFS(mgt,mfermion,mscalar,cc1,cc2)
                    else: Gamma[lab] += gmult[lab]*decayToFS(mgt,mfermion,mscalar,cc1,cc2)
                    
#gravitino -> fermion + vector boson decays
    for lab,pids in fermions.items():
        if not lab in Gamma: Gamma[lab] = 0.
        if lab == 'gluino':
            Gamma[lab] += gmult[lab]*decayToFV(mgt,Masses[pids[0]],0.,2.,0.,0.,0.,0.,0.)
            continue  
        for pid in pids:
            mfermion = Masses[pid]
            iF = (fermions['neutralinos'] + fermions['charginos']).index(pid)        
            for iV,mV in enumerate(mVs[lab]):
                cgg1 = cgl[iV][iF]*conjg(cgl[iV][iF]) + cgr[iV][iF]*conjg(cgr[iV][iF])
                cgg2 = cgl[iV][iF]*conjg(cgr[iV][iF]) + cgr[iV][iF]*conjg(cgl[iV][iF])
                cgh1 = 2.*cgl[iV][iF]*conjg(chL[iV][iF]) - cgr[iV][iF]*conjg(chR[iV][iF])
                cgh2 = 2.*cgl[iV][iF]*conjg(chR[iV][iF]) - cgr[iV][iF]*conjg(chL[iV][iF])
                chh1 = chL[iV][iF]*conjg(chL[iV][iF]) + chR[iV][iF]*conjg(chR[iV][iF])
                chh2 = chL[iV][iF]*conjg(chR[iV][iF]) + chR[iV][iF]*conjg(chL[iV][iF])                
                if mV == mass_W: Gamma[lab] += 2*gmult[lab]*decayToFV(mgt,mfermion,mV,cgg1,cgg2,cgh1,cgh2,chh1,chh2)
                else: Gamma[lab] += gmult[lab]*decayToFV(mgt,mfermion,mV,cgg1,cgg2,cgh1,cgh2,chh1,chh2)         
        
#...g-> axino + axion, axino + saxion:
    if abs(mass_axino) < mgt: Gamma['axino-axion'] = decayToFS(mgt,abs(mass_axino),0.,2.,0.)
    else: Gamma['axino-axion'] = 0.  
    if mass_saxion + abs(mass_axino) < mgt: Gamma['axino-saxion'] = decayToFS(mgt,abs(mass_axino),mass_saxion,2.,0.)
    else: Gamma['axino-saxion'] = 0.    

#Check results
    for gamma in Gamma:
        if abs(Gamma[gamma]) > 0. and Gamma[gamma].imag/abs(Gamma[gamma]) > 10.**-3:
            logger.error("Imaginary width for gravitino -> "+gamma)
            return False
    
#...Get total width and BRs:
    width = abs(sum(Gamma.values())) 
    decays = DecayList()
    decays.width = width
    decays.Xfraction = 1.   
    if width > 0.:
        decays.addDecay(Decay('gravitino',['axino','axion'],abs(Gamma['axino-axion'])/width))
        decays.addDecay(Decay('gravitino',['axino','saxion'],abs(Gamma['axino-saxion'])/width))
        decays.addDecay(Decay('gravitino',['neutralino','radiation'],sum([abs(Gamma[lab]) for lab in labs])/width))
    
    return decays
    
    
    
def decayToFV(mgt,mf,mv,cgg1,cgg2,cgh1,cgh2,chh1,chh2):
    """Master function for computing the decay width for g(mgt) -> fermion(mf) + gauge boson(mv), \
    where fermion = any gaugino/higgsino mixture.
    :param cgg1: real coefficient for the gaugino-gaugino(ll+rr) term
    :param cgg2: real coefficient for the gaugino-gaugino(lr+rl) term
    :param cgh1: real coefficient for the gaugino-higgsino(ll+rr) term
    :param cgh2: real coefficient for the gaugino-higgsino(lr+rl) term
    :param chh1: real coefficient for the higgsino-higgsino(ll+rr) term
    :param chh2: real coefficient for the higgsino-higgsino(lr+rl) term
    """

    if mgt < mv+mf: return 0.
    pq = (mgt**2 + mv**2 - mf**2)/2.    
    pql = (mgt**2 - mv**2 + mf**2)/2.
    qql = (mgt**2 - mv**2 - mf**2)/2.
    msqr = (2./3.)*cgg1*(pq**2*pql/mgt**2 + pq*qql - mv**2*pql)
    msqr += -cgg2*mgt*mf*mv**2
    msqr += (2./3.)*cgh1*mgt*((1./2.)*qql + pq*pql/mgt**2)
    msqr += cgh2*mf*pq
    if mv != 0.:
        msqr += (2./3.)*chh1*(1.+pq**2/(2.*mgt**2*mv**2))*pql
        msqr += (2./3.)*chh2*(1.+pq**2/(2.*mgt**2*mv**2))*mgt*mf
  
    bf = (1./mgt**2)*sqrt(mgt**4 - 2.*mgt**2*(mv**2 + mf**2) + (mf**2 - mv**2)**2)
    return bf*msqr/(16.*Pi*mgt*MPL**2)


def decayToFS(mgt,mf,ms,cc1,cc2):
    """Master function for computing the decay width for g(mgt) -> fermion(mf) + scalar(ms), \
    where fermion = any gaugino/higgsino mixture or lepton,quark and scalar = higgs bosons or slepton,squarks
    :param cc1: real coefficient for the ll+rr term
    :param cc2: real coefficient for the lr+rl term    
    """
        
    if mgt < ms+mf: return 0.
    pq = (mgt**2 + ms**2 - mf**2)/2.
    pql = (mgt**2 - ms**2 + mf**2)/2.     
    msqr = (1./3.)*cc1*(pq**2/mgt**2 - ms**2)*pql
    msqr += (1./3.)*cc2*(pq**2/mgt**2 - ms**2)*mgt*mf
    bf = (1./mgt**2)*sqrt(mgt**4 - 2.*mgt**2*(ms**2 + mf**2) + (mf**2 - ms**2)**2)

    
    return bf*msqr/(16.*Pi*mgt*MPL**2)

def conjg(a):
    """Simple function to return the conjugate value of a complex number."""
     
    if type(a) != type(1j): return a
    else: return a.conjugate()
    
    