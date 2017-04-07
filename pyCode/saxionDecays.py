#!/usr/bin/env python

"""

.. module:: saxionDecays
    :synopsis: Computes saxion BRs and width 

:synopsis: Computes the axino decays
:author: Andre Lessa <lessa.a.p@gmail.com>

"""

import modelParameters
from AuxFuncs import PSlamb, memoize, getParsAt
from AuxDecays import Decay, DecayList
from parameters import Pi,mass_Z,mass_W,sw2,vEW,mt,mb,cw2,mtau
from math import atan,sqrt,atan2,sin,cos,tan

@memoize
def getSaxionBRs(T):
    """Compute axino BRs at temperature T, given the parameters defined in modelParameters.
    Assumes all saxion decay energy goes to radiation except for the decays to axions, axinos or direct decay
    to a pair of lightest neutralinos."""
    
    mass_axino = modelParameters.mass_axino
    mass_saxion = modelParameters.mass_saxion
    Masses = modelParameters.Masses
    slhadata = modelParameters.slhadata
    modelType = modelParameters.modelType
    cH = modelParameters.cH
    vPQ = modelParameters.vPQ
    caYY = modelParameters.caYY
    fa = modelParameters.fa
    nW = modelParameters.nW
    xi = modelParameters.xi
    
    
    

    BRs = {}
    saxGG = saxGLGL = saxAA = saxAT = saxHH = saxZZ = saxWW = saxZH = saxWH = saxFF = saxNN = saxCC = saxSS = 0.
    saxGaGa = saxGaZ = saxN1N1 = 0.    

#Get relevant parameters
    msaxion = abs(mass_saxion)    
    maxino = abs(mass_axino)
    pars = getParsAt(T)
    gSU3 = pars['gSU3']
    gSU2 = pars['gSU2']
    gU1 = pars['gPr']/sqrt(3./5.)
    
    alphaS = gSU3**2/(4.*Pi)
    alphaE = gSU2**2*sw2/(4.*Pi)
    alphaY = alphaE/(1.-sw2)
    alpha2 = alphaE/sw2
    mu = pars['mu']    
    tanb = pars['tanb']
    beta = atan(tanb)  
    mhl = Masses[25]
    mhh = Masses[35]
    ma = Masses[36]
    mhp = Masses[37]


    pidsN = [1000022,1000023,1000025,1000035]
    pidsC = [1000024,1000037]
    mN = [Masses[pid] for pid in pidsN]
    mC = [Masses[pid] for pid in pidsC]
    thZ = [-slhadata["MASS"][pidsN[i]]/mN[i] for i in range(4)]
    thW = [-slhadata["MASS"][pidsC[i]]/mC[i] for i in range(2)]
    v1 = [-slhadata["NMIX"][i,4] for i in range(1,5)]
    v2 = [-slhadata["NMIX"][i,3] for i in range(1,5)]
    v3 = [slhadata["NMIX"][i,2] for i in range(1,5)]
    v4 = [slhadata["NMIX"][i,1] for i in range(1,5)]
    alpha = -slhadata["ALPHA"].values()[0]
    gammaL = atan2(-slhadata["UMIX"][1,1],-slhadata["UMIX"][1,2])
    gammaR = atan2(-slhadata["VMIX"][1,1],-slhadata["VMIX"][1,2])
    thx = slhadata["UMIX"][2,1]/slhadata["UMIX"][1,2]
    thy = slhadata["VMIX"][2,1]/slhadata["VMIX"][1,2]
    thetaT = atan2(slhadata["STOPMIX"][1,2],slhadata["STOPMIX"][1,1])
    thetaB = atan2(slhadata["SBOTMIX"][1,2],slhadata["SBOTMIX"][1,1])
    thetaL = atan2(slhadata["STAUMIX"][1,2],slhadata["STAUMIX"][1,1])
    
    
#...model independent decays:
#..saxion -> gluon + gluon:
    saxGG = alphaS**2*msaxion**3/(32.*Pi**3*fa**2)
#..saxion -> gluino + gluino:
    m = Masses[1000021]
    if msaxion > 2*m: saxGLGL = (alphaS**2*msaxion*m**2/(8.*Pi**3*fa**2))*(1. - 4.*m**2/msaxion**2)**(3./2.)
#...saxion -> axion + axion:
    saxAA = xi**2*msaxion**3/(32.*Pi*fa**2)
#...saxion -> axino + axino:
    if msaxion > 2*maxino: saxAT = xi**2*msaxion*maxino**2 *((1.-4.*maxino**2/msaxion**2)**(3./2.))/(4.*Pi*fa**2)


#...model dependent decays:
    if (modelType == 'DFSZ'):

        lamb = cH*mu/vPQ
        tanr = lamb*vEW*sin(2.*beta)/(2.*mu)
        ephl = -(tanr/(mhl**2 - msaxion**2))*(-ma**2*cos(beta-alpha) + 4.*mu**2*sin(alpha+beta)/sin(2.*beta))
        ephh = -(tanr/(mhh**2 - msaxion**2))*(-ma**2*sin(beta-alpha) + 4.*mu**2*cos(alpha+beta)/sin(2.*beta))
        lshl = (sqrt(2.)*cH*mu**2/vPQ)*(1. - (1./4.)*(ma**2/mu**2)*sin(2.*beta)*sin(2.*alpha)) + (mass_Z**2/(sqrt(2.)*vEW))*cos(2.*alpha)*sin(beta-alpha)*(3.*ephl - ephh*(2.*tan(2.*alpha) + cot(beta-alpha)))                
        lshh = (sqrt(2.)*cH*mu**2/vPQ)*(1. + (1./4.)*(ma**2/mu**2)*sin(2.*beta)*sin(2.*alpha)) + (mass_Z**2/(sqrt(2.)*vEW))*cos(2.*alpha)*sin(beta-alpha)*(3.*ephh*cot(beta-alpha) + ephl*(2.*tan(2.*alpha)*cot(beta-alpha) -1.))
        lshlhh = -(cH*ma**2/(2.*sqrt(2.)*vPQ))*sin(2.*beta)*cos(2.*alpha) + (mass_Z**2/(sqrt(2.)*vEW))*cos(2.*alpha)*sin(beta-alpha)*(-ephl*(2.*tan(2.*alpha) + cot(beta-alpha))+ ephh*(2.*tan(2.*alpha)*cot(beta-alpha) - 1.))
        lsa = (sqrt(2.)*cH*mu**2/vPQ)*(1. + (1./4.)*(ma**2/mu**2)*sin(2.*beta)**2) + (mass_Z**2/(sqrt(2.)*vEW))*cos(2.*beta)*sin(beta-alpha)*(ephl - ephh*cot(beta-alpha))
        lshphm = (sqrt(2.)*cH*mu**2/vPQ)*(1. + (1./4.)*(ma**2/mu**2)*sin(2.*beta)**2) + (mass_Z**2/(sqrt(2.)*vEW))*(cos(2.*beta)*sin(beta-alpha)*(ephl - ephh*cot(beta-alpha))+ 2.*cw2*sin(beta+alpha)*(ephl + ephh*cot(beta+alpha)))
        gsvv = ephl*sin(beta+alpha) + ephh*cos(beta+alpha)
        gsvh = ephl*cos(beta+alpha) - ephh*sin(beta+alpha)
        ggz = gSU2/sqrt(cw2)
        ggw = gSU2
        gsffu = (1./sin(beta))*(-ephl*cos(alpha) + ephh*sin(alpha))
        gsffd = (1./cos(beta))*(-ephl*sin(alpha) - ephh*cos(alpha))
        xhl,xhh,xa,xs = [0.]*len(mN),[0.]*len(mN),[0.]*len(mN),[0.]*len(mN)
        for iN,thZi in enumerate(thZ):
            xhl[iN],xhh[iN],xa[iN],xs[iN] = [0.]*len(mN),[0.]*len(mN),[0.]*len(mN),[0.]*len(mN)
            for jN,thZj in enumerate(thZ):
                xhl[iN][jN] = -(1./2.)*thZi*thZj*(v2[iN]*sin(alpha) - v1[iN]*cos(alpha))*(gSU2*v3[jN] - gU1*v4[jN])
                xhh[iN][jN] = -(1./2.)*thZi*thZj*(v2[iN]*cos(alpha) + v1[iN]*sin(alpha))*(gSU2*v3[jN] - gU1*v4[jN])
                xa[iN][jN] = (1./2.)*thZi*thZj*(v2[iN]*sin(beta) - v1[iN]*cos(beta))*(gSU2*v3[jN] - gU1*v4[jN])
                xs[iN][jN] = thZi*thZj*v1[iN]*v2[jN]
        s1hl = (1./2.)*thW[0]*(sin(alpha)*sin(gammaR)*cos(gammaL) + cos(alpha)*sin(gammaL)*cos(gammaR))
        s2hl = -(1./2.)*thW[1]*thx*thy*(sin(alpha)*cos(gammaR)*sin(gammaL) + cos(alpha)*cos(gammaL)*sin(gammaR))
        shl = (1./2.)*(-thW[0]*thx*sin(alpha)*sin(gammaR)*sin(gammaL) + thW[0]*thx*cos(alpha)*cos(gammaL)*cos(gammaR)- thW[1]*thy*cos(alpha)*sin(gammaR)*sin(gammaL) + thW[1]*thy*sin(alpha)*cos(gammaR)*cos(gammaL))
        phl = (1./2.)*(thW[0]*thx*sin(alpha)*sin(gammaR)*sin(gammaL) - thW[0]*thx*cos(alpha)*cos(gammaL)*cos(gammaR)- thW[1]*thy*cos(alpha)*sin(gammaR)*sin(gammaL) + thW[1]*thy*sin(alpha)*cos(gammaR)*cos(gammaL))
        s1hh = (1./2.)*thW[0]*(cos(alpha)*sin(gammaR)*cos(gammaL) - sin(alpha)*sin(gammaL)*cos(gammaR))
        s2hh = -(1./2.)*thW[0]*thx*thy*(cos(alpha)*cos(gammaR)*sin(gammaL) - sin(alpha)*cos(gammaL)*sin(gammaR))
        shh = (1./2.)*(-thW[0]*thx*cos(alpha)*sin(gammaR)*sin(gammaL) - thW[0]*thx*sin(alpha)*cos(gammaL)*cos(gammaR)+ thW[1]*thy*sin(alpha)*sin(gammaR)*sin(gammaL) + thW[1]*thy*cos(alpha)*cos(gammaR)*cos(gammaL))
        phh = (1./2.)*(thW[0]*thx*cos(alpha)*sin(gammaR)*sin(gammaL) + thW[0]*thx*sin(alpha)*cos(gammaL)*cos(gammaR)+ thW[1]*thy*sin(alpha)*sin(gammaR)*sin(gammaL) + thW[1]*thy*cos(alpha)*cos(gammaR)*cos(gammaL))
        s1s = thW[0]*cos(gammaL)*cos(gammaR)
        s2s = -thW[1]*thx*thy*sin(gammaL)*sin(gammaR)
        ss = (1./2.)*(thW[0]*thy*cos(gammaL)*sin(gammaR)- thW[1]*thx*sin(gammaL)*cos(gammaR))
        ps = (1./2.)*(thW[0]*thy*cos(gammaL)*sin(gammaR)+ thW[1]*thx*sin(gammaL)*cos(gammaR))
        sig1s = ephl*s1hl + ephh*s1hh - cH*mu*s1s/(4.*gSU2*vPQ)
        sig2s = ephl*s2hl + ephh*s2hh - cH*mu*s2s/(4.*gSU2*vPQ)
        sigs = ephl*shl + ephh*shh - cH*mu*ss/(2.*gSU2*vPQ)
        Pis = ephl*phl + ephh*phh - cH*mu*ps/(2.*gSU2*vPQ)
        zets = [0.]*len(mN)
        for iN,xsi in enumerate(xs):
            zets[iN] = [ephl*xhl[iN][jN] + ephh*xhh[iN][jN] - cH*mu*xsij/(2.*sqrt(2.)*vPQ) for jN,xsij in enumerate(xsi)]
        
        labSQUP = ['sq_up','sq_ch','sq_tp']
        mSQUP = [[Masses[1000002],Masses[2000002]],[Masses[1000004],Masses[2000004]],[Masses[1000006],Masses[2000006]]]
        mqUP = [0.,0.,mt]
        labSQDW = ['sq_dn','sq_st','sq_bt']
        mSQDW = [[Masses[1000001],Masses[2000001]],[Masses[1000003],Masses[2000003]],[Masses[1000005],Masses[2000005]]]
        mqDW = [0.,0.,mb]
        labSL = ['sl_e','sl_m','sl_t']
        mSL = [[Masses[1000011],Masses[2000011]],[Masses[1000013],Masses[2000013]],[Masses[1000015],Masses[2000015]]]
        ml = [0.,0.,mtau]
        labSN = ['nue','num','nut']
        mSN = [[Masses[1000012]],[Masses[1000014]],[Masses[1000016]]]
        mnu = [0.,0.,0.]
        aUP = slhadata["AU"].values()
        aDW = slhadata["AD"].values()
        aSL = slhadata["AE"].values()
        aSN = [0.]*len(mSN)
        thetaUP = [0.,0.,thetaT]
        thetaDW = [0.,0.,thetaB]
        thetaSL = [0.,0.,thetaL]
        thetaSN = [0.]*len(mSN)
        labSF = labSQUP + labSQDW + labSL + labSN
        mSF = mSQUP + mSQDW + mSL + mSN
        mF = mqUP + mqDW + ml + mnu
        aSF = aUP + aDW + aSL + aSN
        thetaSF = thetaUP + thetaDW + thetaSL + thetaSN
        FHcouplings = {'H': {}, 'h' : {}} 
        for iF,labF in enumerate(labSF):   #sfermions                        
            mf = mF[iF]
            A = aSF[iF]
            thetaF = thetaSF[iF]
            cosF = cos(thetaF)
            sinF = sin(thetaF)            
            if labF in labSQUP+labSN:  #up squarks + sneutrinos 
                ahll = gSU2*(mass_W*(1./2. - (1./6.)*(sw2/cw2))*sin(beta-alpha) - (mf**2/mass_W)*cos(alpha)/sin(beta))
                aHll = gSU2*(-mass_W*(1./2. - (1./6.)*(sw2/cw2))*cos(beta-alpha) + (mf**2/mass_W)*sin(alpha)/sin(beta))
                if labF in labSQUP:
                    ahrr = gSU2*((2./3.)*mass_W*(sw2/cw2)*sin(beta-alpha) - (mf**2/mass_W)*cos(alpha)/sin(beta))
                    ahlr = gSU2*(mf/(2.*mass_W*sin(beta)))*(-mu*sin(alpha) + A*cos(alpha))
                    aHrr = gSU2*(-(2./3.)*mass_W*(sw2/cw2)*cos(beta-alpha) + (mf**2/mass_W)*sin(alpha)/sin(beta))
                    aHlr = gSU2*(mf/(2.*mass_W*sin(beta)))*(-mu*cos(alpha) - A*sin(alpha))
                else: ahrr = ahlr = aHrr = aHlr = 0.  #sneutrinos
            elif labF in labSQDW+labSL:  #down squarks + sleptons
                ahll = gSU2*(mass_W*(-1./2. - (1./6.)*(sw2/cw2))*sin(beta-alpha) - (mf**2/mass_W)*sin(alpha)/cos(beta))
                ahrr = gSU2*((-1./3.)*mass_W*(sw2/cw2)*sin(beta-alpha) - (mf**2/mass_W)*sin(alpha)/cos(beta))
                ahlr = gSU2*(mf/(2.*mass_W*cos(beta)))*(-mu*cos(alpha) + A*sin(alpha))
                aHll = gSU2*(mass_W*(1./2. + (1./6.)*(sw2/cw2))*cos(beta-alpha) - (mf**2/mass_W)*cos(alpha)/cos(beta))
                aHrr = gSU2*((1./3.)*mass_W*(sw2/cw2)*cos(beta-alpha) - (mf**2/mass_W)*cos(alpha)/cos(beta))
                aHlr = gSU2*(mf/(2.*mass_W*cos(beta)))*(mu*sin(alpha) + A*cos(alpha))
            aHiggs = {'h' : [ahll,ahlr,ahrr], 'H' : [aHll,aHlr,aHrr]}          
            for labH, in aHiggs:                
                aF = aHiggs[labH][:]
                FHcouplings[labH][labF] = []                    
                FHcouplings[labH][labF].append([aF[0]*cosF**2 + aF[2]*sinF**2 - 2.*aF[1]*cosF*sinF
                              ,aF[0]*cosF*sinF - aF[2]*cosF*sinF + 2.*aF[1]*cos(2.*thetaF)])
                FHcouplings[labH][labF].append([aF[0]*cosF*sinF - aF[2]*cosF*sinF + 2.*aF[1]*cos(2.*thetaF)
                              ,aF[0]*sinF**2 + aF[2]*cosF**2 + 2.*aF[1]*cosF*sinF])

#...dfsz decays:
#...decays to higgs:
        lsall = []
        lsall.append([lshl,lshlhh,0.,0.,0.])
        lsall.append([lshlhh,lshh,0.,0.,0.])
        lsall.append([0.,0.,lsa,0.,0.])
        lsall.append([0.,0.,0.,0.,lshphm])
        lsall.append([0.,0.,0.,lshphm,0.])
        mhv = [mhl,mhh,ma,mhp,mhp]
        for ih,mhi in enumerate(mhv):
            for jh,mhj in enumerate(mhv):
                if msaxion > mhi + mhj:
                    PS = sqrt(PSlamb(1.,mhi**2/msaxion**2,mhj**2/msaxion**2))
                    saxHH += (lsall[ih][jh]**2/(16.*Pi*msaxion))*PS/2.
        
#...decays to vector boson pair:
        if msaxion > 2*mass_Z:
            saxZZ = (ggz**2*gsvv**2/(16.*Pi))*msaxion*(3.*mass_Z**2/msaxion**2 
                     + (msaxion**2/(4.*mass_Z**2)))*(1. - 4.*mass_Z**2/msaxion**2)*sqrt(1. - 4.*mass_Z**2/msaxion**2)/2.
        if msaxion > 2*mass_W:
            saxWW = (ggw**2*gsvv**2/(16.*Pi))*msaxion*(3.*mass_W**2/msaxion**2 
                     + (msaxion**2/(4.*mass_W**2)))*(1. - 4.*mass_W**2/msaxion**2)*sqrt(1. - 4.*mass_W**2/msaxion**2)

#...decays to vector boson - higgs scalar:
        if msaxion > mass_Z + ma:
            PS = sqrt(PSlamb(1.,mass_Z**2/msaxion**2,ma**2/msaxion**2))
            saxZH = (ggz**2*gsvh**2/(32.*Pi))*(msaxion**3/mass_Z**2)*((1. - ma**2/msaxion**2)**2 
                                    - 2.*(mass_Z**2/msaxion**2)*(1. + ma**2/msaxion**2) + mass_Z**4/msaxion**4)*PS/2.
        if msaxion > mass_W + mhp:
            PS = sqrt(PSlamb(1.,mass_W**2/msaxion**2,mhp**2/msaxion**2))
            saxWH = (ggw**2*gsvh**2/(32.*Pi))*(msaxion**3/mass_W**2)*((1. - mhp**2/msaxion**2)**2 
                                    - 2.*(mass_W**2/msaxion**2)*(1. + mhp**2/msaxion**2) + mass_W**4/msaxion**4)*PS/2.

#...decay to SM fermions:
#...top-top:
        if msaxion > 2*mt:
            saxFF += (3./(16.*Pi))*(mt**2/vEW**2)*gsffu**2*msaxion*(1. - 4.*mt**2/msaxion**2)**(3./2.)
        if msaxion > 2*mb:
            saxFF += (3./(16.*Pi))*(mb**2/vEW**2)*gsffd**2 *msaxion*(1. - 4.*mb**2/msaxion**2)**(3./2.)

#...decay to chargino pair:
        allsigs = [sig1s,sig2s]
        for iC,mc in enumerate(mC):
            if msaxion > 2*mc: saxCC += (gSU2**2/(4.*Pi))*allsigs[iC]**2*msaxion*(1. - 4.*mc**2/msaxion**2)**(3./2.)

        if msaxion > sum(mC):
            PS = sqrt(PSlamb(1.,mC[0]**2/msaxion**2,mC[1]**2/msaxion**2))
            saxCC += 2.*(gSU2**2/(16.*Pi))*msaxion*PS*(sigs**2*(1. - (sum(mC)/msaxion)**2) 
                                + Pis**2*(1. - ((mC[1]-mC[0])/msaxion)**2))     # the factor of 2 takes care of w1w2+w2w1

#...decay to neutralino pair:
        for iN,mni in enumerate(mN):
            for jN,mnj in enumerate(mN):
                if msaxion > mni+mnj:
                    PS = sqrt(PSlamb(1.,mni**2/msaxion**2,mnj**2/msaxion**2))
                    saxnn = (1./(8.*Pi))*msaxion*(zets[iN][jN]+zets[jN][iN])**2*(1. 
                                    - ((mni + thZ[iN]*thZ[jN]*mnj)/msaxion)**2)*PS/2.  
                    if iN == jN == 0: saxN1N1 = saxnn
                    else: saxNN += saxnn
                    
#...decay to sfermion + sfermion
        for iF,labF in enumerate(labSF):
            if labF in labSQUP + labSQDW: ncf = 3.
            else: ncf = 1.            
            msfv = mSF[iF]
            for jf,msfj in enumerate(msfv):
                for kf,msfk in enumerate(msfv):
                    if msaxion > msfj+msfk:
                        PS = sqrt(PSlamb(1.,msfj**2/msaxion**2,msfk**2/msaxion**2))
                        ahf = FHcouplings['h'][labF][jf][kf]
                        aHf = FHcouplings['H'][labF][jf][kf]
                        saxSS += (1./(16.*Pi*msaxion))*((ephl*ahf + ephh*aHf)**2*ncf*PS)

#--------------------------------------------------------
#...ksvz decays:
    elif (modelType == 'KSVZ'):
        if nW != 2: alpha2 = 0.

#...saxion -> w+w
        if msaxion > 2.*mass_W:
            gw = (alpha2/(32.*Pi))*3/fa
            PS = sqrt(1. - 4.*mass_W**2/msaxion**2)
            saxWW = 2.*gw**2*msaxion**3*PS*(1.- 4.*mass_W**2/msaxion**2 + 6.*mass_W**4/msaxion**4)/Pi

#...saxion -> z+z
        if msaxion > 2.*mass_Z:
            gz = (alpha2/(32.*Pi))*3*cw2/fa  + (alphaY/(16.*Pi))*caYY*sw2/fa
            PS = sqrt(1. - 4.*mass_Z**2/msaxion**2)
            saxZZ = gz**2*msaxion**3*PS*(1. - 4.*mass_Z**2/msaxion**2 + 6.*mass_Z**4/msaxion**4)/Pi

#...saxion -> gamma+gamma
        gg = (alpha2/(32.*Pi))*3*sw2/fa  + (alphaY/(16.*Pi))*caYY*cw2/fa
        saxGaGa = gg**2*msaxion**3/Pi
#...saxion -> gamma+z
        if msaxion > mass_Z:
            gz = (-(alpha2/(32.*Pi))*3 + (alphaY/(16.*Pi))*caYY)*sqrt(sw2*cw2)/fa
            saxGaZ = 2.*gz**2*msaxion**3*(1.-mass_Z**2/msaxion**2)**4/Pi

#...saxion -> chargino+chargino                                             # obs: i=1,j=2 is equal i=2,j=1 (w_1^- + w_2^+)#
        ql = [[cos(gammaL)**2,abs(cos(gammaL))*sin(gammaL)],[abs(cos(gammaL))*sin(gammaL),sin(gammaL)**2]]
        qr = [[sin(gammaR)**2,abs(cos(gammaR))*sin(gammaR)],[abs(cos(gammaR))*sin(gammaR),cos(gammaR)**2]]
        for iC,mci in enumerate(mC):
            for jC,mcj in enumerate(mC):                           
                if msaxion > mci+mcj:                    
                    gw = (alpha2/(32.*Pi))*3/fa
                    PS = sqrt(PSlamb(1.,mci**2/msaxion**2,mcj**2/msaxion**2)) 
                    saxCC += gw**2*msaxion*PS*((ql[iC][jC]**2 + qr[iC][jC]**2)*((mci**2+mcj**2)*(1. 
                       - mci**2/msaxion**2 - mcj**2/msaxion**2) - 4.*mci**2*mcj**2/msaxion**2) 
                       + (2.*ql[iC][jC]*qr[iC][jC])*2.*mci*mcj*(1. - 2.*mci**2/msaxion**2 - 2.*mcj**2/msaxion**2))/Pi

#...saxion -> neutralino+neutralino:
        for iN,mni in enumerate(mN):
            for jN,mnj in enumerate(mN):
                if msaxion > mni+mnj:
                    gz = (alpha2/(32.*Pi))*3*v3[iN]*v3[jN]/fa  + (alphaY/(16.*Pi))*(caYY)*v4[iN]*v4[jN]/fa
                    PS = sqrt(PSlamb(1.,mni**2/msaxion**2,mnj**2/msaxion**2))
                    saxNN += gz**2*msaxion*PS*((mni + mnj)**2)*(1. - (mni + mnj)**2/msaxion**2)/Pi


#...Get total width and BRs:
#   axodd = axnh + axch + tgl
    Gamma = {} 
    Gamma['saxion-Gl-Gl'] = saxGLGL
    Gamma['saxion-g-g'] = saxGG
    Gamma['saxion-axion-axion'] = saxAA
    Gamma['saxion-axino-axino'] = saxAT
    Gamma['saxion-H-H'] = saxHH
    Gamma['saxion-gamma-gamma'] = saxGaGa
    Gamma['saxion-gamma-Z'] = saxGaZ
    Gamma['saxion-Z-Z'] = saxZZ
    Gamma['saxion-W-W'] = saxWW
    Gamma['saxion-Z-H'] = saxZH
    Gamma['saxion-W-H'] = saxWH
    Gamma['saxion-f-f'] = saxFF
    Gamma['saxion-N-N'] = saxNN
    Gamma['saxion-N1-N1'] = saxN1N1
    Gamma['saxion-C-C'] = saxCC
    Gamma['saxion-sf-sf'] = saxSS        
    width = sum(Gamma.values())
    decays = DecayList()
    decays.width = width

        
    if width > 0.:
        decays.Xfraction = 1.-(saxAA+saxAT)/width
        decays.addDecay(Decay('saxion',['radiation','radiation'],
                              (saxGG+saxHH+saxGaGa+saxGaZ+saxZZ+saxWW+saxZH+saxWH+saxFF)/width))
        decays.addDecay(Decay('saxion',['neutralino','neutralino'],(saxNN+saxN1N1+saxCC+saxSS+saxGLGL)/width))
        decays.addDecay(Decay('saxion',['axion','axion'],saxAA/width))
        decays.addDecay(Decay('saxion',['axino','axino'],saxAT/width))

    return decays


def cot(x):
    """Cotangent function"""
    
    return 1./tan(x)
