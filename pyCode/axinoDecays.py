#!/usr/bin/env python

"""

.. module:: axinoDecays
    :synopsis: Computes axino BRs and width 

:synopsis: Computes the axino decays
:author: Andre Lessa <lessa.a.p@gmail.com>, Hasan Baris

"""

import modelParameters
from AuxFuncs import PSlamb, memoize, getParsAt
from AuxDecays import Decay, DecayList
from parameters import Pi,mass_Z,mass_W,sw2,vEW,mt,mb
from math import atan,sqrt,atan2,sin,cos

# @memoize
def getAxinoBRs(T):
    """Compute axino BRs at temperature T, given the parameters defined in modelParameters.\
    Assumes all axino decay energy goes into radiation."""
    
    mass_axino = modelParameters.mass_axino
    Masses = modelParameters.Masses
    slhadata = modelParameters.slhadata
    modelType = modelParameters.modelType
    cH = modelParameters.cH
    vPQ = modelParameters.vPQ
    caYY = modelParameters.caYY
    fa = modelParameters.fa
    nW = modelParameters.nW
    
    axGL = axCH = axNZ = axNG = axNH = axCW = axFSF = 0.
    
#Get relevant parameters  
    maxino = abs(mass_axino)    
    pars = getParsAt(T)
    gSU3 = pars['gSU3']
    gSU2 = pars['gSU2']
    gU1 = pars['gPr']
    alphaS = gSU3**2/(4.*Pi)
    alphaE = gSU2**2*sw2/(4.*Pi)
    alphaY = alphaE/(1.-sw2)
    alpha2 = alphaE/sw2
    mu = pars['mu']    
    tanb = pars['tanb']
    beta = atan(tanb)
    sq2=float(sqrt(2.))    
    mhl = Masses[25]
    mhh = Masses[35]
    ma = Masses[36]
    mhp = Masses[37]
    stmass = [Masses[1000006],Masses[2000006]]
    sbmass = [Masses[1000005],Masses[2000005]]    

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
    Vs = [v1,v2,v3,v4]
    alpha = -slhadata["ALPHA"].values()[0]
    gammaL = atan2(-slhadata["UMIX"][1,1],-slhadata["UMIX"][1,2])
    gammaR = atan2(-slhadata["VMIX"][1,1],-slhadata["VMIX"][1,2])
    thx = slhadata["UMIX"][2,1]/slhadata["UMIX"][1,2]
    thy = slhadata["VMIX"][2,1]/slhadata["VMIX"][1,2]
    thetaT = atan2(slhadata["STOPMIX"][1,2],slhadata["STOPMIX"][1,1])
    thetaB = atan2(slhadata["SBOTMIX"][1,2],slhadata["SBOTMIX"][1,1])
    
    if mass_axino > 0.:
        thaxi = 1.
        thax = 1.
    else:
        thaxi = 1.j
        thax = -1.
        

#...model independent decays:
#...axino-> gluino + gluon    
    m = Masses[1000021]    
    if maxino > m: axGL = 8.*alphaS**2*maxino**3*(1.-(m/maxino)**2)**3/(128.*Pi**3*fa**2)
    else: axGL = 0. 
    
 
#...model dependent decays:
    if modelType == 'DFSZ':

        tazhl = [thZ[i]*thax*(v1[i]*sin(alpha) + v2[i]*cos(alpha)) for i in range(4)]
        tazhh = [thZ[i]*thax*(v1[i]*cos(alpha) - v2[i]*sin(alpha)) for i in range(4)]        
        taza = [-thZ[i]*thax*(-v1[i]*sin(beta) + v2[i]*cos(beta)) for i in range(4)] 
        singr = sin(gammaR)
        cosgr = cos(gammaR)
        singl = sin(gammaL)
        cosgl = cos(gammaL)
        v0 = [sum([(cH*mu*vEW/vPQ)*Vs[iv][i]*(v1[i]*cos(beta)+v2[i]*sin(beta))/(thax*maxino-mN[i]*thZ[i]) for i in range(4)]) for iv in range(4)]
 
#     neutralino couplings
        xchl1 = [(-1./2.)*thZ[i]*thax*(v2[i]*sin(alpha)-v1[i]*cos(alpha))*(gSU2*v0[2]-gU1*v0[3]) for i in range(4)]
        xchl2 = [(-1./2.)*thZ[i]*thax*(v0[1]*sin(alpha)-v0[0]*cos(alpha))*(gSU2*v3[i]-gU1*v4[i]) for i in range(4)]
        xchh1 = [(-1./2.)*thZ[i]*thax*(v2[i]*cos(alpha)+v1[i]*sin(alpha))*(gSU2*v0[2]-gU1*v0[3]) for i in range(4)]
        xchh2 = [(-1./2.)*thZ[i]*thax*(v0[1]*cos(alpha)+v0[0]*sin(alpha))*(gSU2*v3[i]-gU1*v4[i]) for i in range(4)]
        xca1 = [(1./2.)*thZ[i]*thax*(v2[i]*sin(beta)-v1[i]*cos(beta))*(gSU2*v0[2]-gU1*v0[3]) for i in range(4)]
        xca2 = [(1./2.)*thZ[i]*thax*(v0[1]*sin(beta)-v0[0]*cos(beta))*(gSU2*v3[i]-gU1*v4[i]) for i in range(4)]

        lmbazhl = [xchl1[i]+xchl2[i]-cH*mu*tazhl[i]/(sq2*vPQ) for i in range(4)]     
        lmbazhh = [xchh1[i]+xchh2[i]-cH*mu*tazhh[i]/(sq2*vPQ) for i in range(4)]
        lmbaza = [xca1[i]+xca2[i]-cH*mu*taza[i]/(sq2*vPQ) for i in range(4)]
         
        a0 = [0.]*len(mN)
        a0[0] = -(1./sq2)*(gSU2*v0[2]+gU1*v0[3])*singr-gSU2*v0[0]*cosgr
        a0[1] = (1./sq2)*(gSU2*v0[2]+gU1*v0[3])*cosgr-gSU2*v0[0]*singr
        a0[2] = -(1./sq2)*(gSU2*v0[2]+gU1*v0[3])*singl+gSU2*v0[1]*cosgl
        a0[3] = (1./sq2)*(gSU2*v0[2]+gU1*v0[3])*cosgl+gSU2*v0[1]*singl
        lmb10 = a0[0]+thax*cH*mu*tanb*singr/vPQ
        lmb20 = a0[1]-thax*cH*mu*tanb*cosgr/vPQ
        lmb30 = a0[2]-thax*cH*mu/tanb*singl/vPQ
        lmb40 = a0[3]+thax*cH*mu/tanb*cosgl/vPQ
        a1 = (1./2.)*(thW[0]*cos(beta)*lmb20 - thax*sin(beta)*lmb40)
        b1 = (1./2.)*(thW[0]*cos(beta)*lmb20 + thax*sin(beta)*lmb40)
        a2 = (1./2.)*(thW[1]*thy*cos(beta)*lmb10-thax*thx*sin(beta)*lmb30)
        b2 = (1./2.)*(thW[1]*thy*cos(beta)*lmb10+thax*thx*sin(beta)*lmb30)
        x0, y0 = [0.]*len(mC),[0.]*len(mC)
        x0[0] = (1./2.)*(thW[0]*thax*(cosgr*v0[0]/sq2+singr*v0[2])-cosgl*v0[1]/sq2+singl*v0[2])
        x0[1] = (1./2.)*(thW[1]*thax*thy*(-singr*v0[0]/sq2+cosgr*v0[2])+thx*(singl*v0[1]/sq2+cosgl*v0[2]))
        y0[0] = (1./2.)*(-thW[0]*thax*(cosgr*v0[0]/sq2+singr*v0[2])-cosgl*v0[1]/sq2+singl*v0[2])
        y0[1] = (1./2.)*(-thW[1]*thax*thy*(-singr*v0[0]/sq2+cosgr*v0[2])+thx*(singl*v0[1]/sq2+cosgl*v0[2]))
        wio = [0.]*len(mN)
        for i in range(4):            
            wio[i] = (1./4.)*sqrt(gSU2**2+gU1**2)*(v1[i]*v0[0]-v2[i]*v0[1])
            if thZ[i] > 0.: wio[i] *= thaxi
 
        ytop = mt*sq2/vEW
        ybot = mb*sq2/vEW
        axu = (-1.*(-thaxi)/(-1j*sq2))*(gSU2*v0[2]+gU1/3.*v0[3])
        axd = (-1.*(-thaxi)/(-1j*sq2))*(-gSU2*v0[2]+gU1/3.*v0[3])
        bxu = (4./(3.*sq2)*thaxi/1j)*gU1*v0[3]
        bxd = (-2./(3.*sq2)*thaxi/1j)*gU1*v0[3]        
        au = [0.]*len(mC)
        bu = [0.]*len(mC)
        ad = [0.]*len(mC)
        bd = [0.]*len(mC)
        au[0] = (1./2.)*((1j*axu-thaxi*ytop*v0[0])*cos(thetaT)-(1j*bxu-(-thaxi)*ytop*v0[0])*sin(thetaT))
        bu[0] = (1./2.)*((-1j*axu-thaxi*ytop*v0[0])*cos(thetaT)-(1j*bxu+(-thaxi)*ytop*v0[0])*sin(thetaT))
        ad[0] = (1./2.)*((1j*axd-thaxi*ybot*v0[1])*cos(thetaB)-(1j*bxd-(-thaxi)*ybot*v0[1])*sin(thetaB))
        bd[0] = (1./2.)*((-1j*axd-thaxi*ybot*v0[1])*cos(thetaB)-(1j*bxd+(-thaxi)*ybot*v0[1])*sin(thetaB))
        au[1] = (1./2.)*((1j*axu-thaxi*ytop*v0[0])*sin(thetaT)+(1j*bxu-(-thaxi)*ytop*v0[0])*cos(thetaT))
        bu[1] = (1./2.)*((-1j*axu-thaxi*ytop*v0[0])*sin(thetaT)+(1j*bxu+(-thaxi)*ytop*v0[0])*cos(thetaT))
        ad[1] = (1./2.)*((1j*axd-thaxi*ybot*v0[1])*sin(thetaB)+(1j*bxd-(-thaxi)*ybot*v0[1])*cos(thetaB))
        bd[1] = (1./2.)*((-1j*axd-thaxi*ybot*v0[1])*sin(thetaB)+(1j*bxd+(-thaxi)*ybot*v0[1])*cos(thetaB))


#...dfsz decays:
#...decays to neutralino + higgs:
        mHs = [mhl,mhh,ma]
        lmbazs = [lmbazhl,lmbazhh,lmbaza]
        hsign = [1,1,-1]
        for ih,mh in enumerate(mHs):            
            for iN,mn in enumerate(mN):
                if maxino > mn+mh:                                  
                    PS = sqrt(PSlamb(1.,mn**2/maxino**2,mh**2/maxino**2))
                    axNH += (1./(16.*Pi))*lmbazs[ih][iN]**2*maxino*PS*((1. + mn**2/maxino**2 - mh**2/maxino**2) 
                                                                   + hsign[ih]*2.*thZ[iN]*thax*mn/maxino)
 
#...decays to chargino + higgs:
        a = [a1,a2]
        b = [b1,b2]
        for iC,mc in enumerate(mC):
            if maxino > mc+mhp:
                PS = sqrt(PSlamb(1.,mc**2/maxino**2,mhp**2/maxino**2))
                axCH += 2.*(1./(16.*Pi))*maxino*PS*((a[iC]**2 + b[iC]**2)*(1. + mc**2/maxino**2 - mhp**2/maxino**2) 
                                                + 2.*(a[iC]**2 - b[iC]**2)*mc/maxino)     	#the factor of two takes care of both charge combinations
                  
 
#...decays to neutralino + z-boson:                               -- from arxiv:1309.5365,hasan
        for iN,mn in enumerate(mN):
            if maxino > mn+mass_Z:
                PS = sqrt(PSlamb(1.,mn**2/maxino**2,mass_Z**2/maxino**2))
                axNZ += (1./(4.*Pi))*abs(wio[iN]**2)*maxino*PS*((1. + mn**2/maxino**2 - 2*(mass_Z**2/maxino**2)) 
                                                            + (maxino**2/mass_Z**2)*(1. - mn**2/maxino**2)**2 + 6*thZ[iN]*thax*mn/maxino)
                   
        
#...decays to chargino + w-boson:                                 -- from arxiv:1309.5365,hasan
        for iC,mc in enumerate(mC):        
            if maxino > mc+mass_W:
                PS = sqrt(PSlamb(1.,mc**2/maxino**2,mass_W**2/maxino**2))
                axCW += 2.*(1./(16.*Pi))*maxino*gSU2**2*PS*((x0[iC]**2+y0[iC]**2)*((1+mc**2/maxino**2-2*mass_W**2/maxino**2)
                                        +(maxino**2/mass_W**2)*(1-mc**2/maxino**2)**2)-6*(x0[iC]**2-y0[iC]**2)*mc/maxino)   #the factor of two takes care of both charge combinations
 
        
#...decays to fermion + sfermion:                                 -- from arxiv:1309.5365,hasan    
        for iF,msf in enumerate(stmass):
            if maxino > msf+mt:
                PS = sqrt(PSlamb(1.,msf**2/maxino**2,mt**2/maxino**2))
                axFSF += 3.*maxino/(16.*Pi)*PS*(abs(au[iF]**2)*((1.+mt/maxino)**2
                                            -(msf**2/maxino**2))+abs(bu[iF]**2)*((1.-mt/maxino)**2-(msf**2/maxino**2)))
        for iF,msf in enumerate(sbmass):
            if maxino > msf+mb:
                PS = sqrt(PSlamb(1.,msf**2/maxino**2,mb**2/maxino**2))
                axFSF += 3.*maxino/(16.*Pi)*PS*(abs(ad[iF]**2)*((1.+mb/maxino)**2
                                            -(msf**2/maxino**2))+abs(bd[iF]**2)*((1.-mb/maxino)**2-(msf**2/maxino**2)))

#-------------------------------------------------------------------------
#...ksvz decays:
    elif modelType == 'KSVZ':
#...axino -> neutralino_i + gamma,Z
        for iN,mn in enumerate(mN):
            gg = (alphaY/(16.*Pi))*(caYY)*v4[iN]*sqrt(1.-sw2)/fa  # coupling for singlet pq fermions
            gz = (alphaY/(16.*Pi))*(caYY)*v4[iN]*sqrt(sw2)/fa
            if nW == 2 :
                gg += (alpha2/(32.*Pi))*3*v3[iN]*sqrt(sw2)/fa  # coupling for doublet pq fermions      
                gz += - (alpha2/(32.*Pi))*3*v3[iN]*sqrt(1.-sw2)/fa       
            if mn < maxino: axNG += 2.*gg**2*((maxino**2 - mn**2)**3)/(Pi*maxino**3)
            if mn+mass_Z < maxino:
                PS = sqrt(PSlamb(1.,mn**2/maxino**2,mass_Z**2/maxino**2))
                axNZ += 2.*gz**2*maxino**3 *PS*((1. - mn**2/maxino**2)**2 + 3.*mn*mass_Z**2/maxino**3 
                                        - (mass_Z**2/(2.*maxino**2))*(1. + mn**2/maxino**2 + mass_Z**2/maxino**2))/Pi
 
#...axino -> chargino_i + w
        gammaM = [[sin(gammaL),sin(gammaR)],[abs(cos(gammaL)),abs(cos(gammaR))]]
        for iC,mc in enumerate(mC):
            if nW != 2: continue      # only for doublet pq fermions :
            uu,vv = gammaM[iC]
            gw = (alpha2/(32.*Pi))*3/fa
            if mc+mass_W < maxino:
                PS = sqrt(PSlamb(1.,mc**2/maxino**2,mass_W**2/maxino**2))
                axCW += 2.*gw**2*maxino**3*PS*(((1. - mc**2/maxino**2)**2 
                                        - (mass_W**2/(2.*maxino**2))*(1. + mc**2/maxino**2 
                                        + mass_W**2/maxino**2))*(uu**2+vv**2) + 3.*(mc*mass_W**2/maxino**3)*(2.*uu*vv))/Pi
 
 
 
#...Get total width and BRs:
#   axodd = axnh + axch + tgl
    Gamma = {}
    Gamma['axino-Gl-gl'] = axGL
    Gamma['axino-N-higgs'] = axNH
    Gamma['axino-C-higgs'] = axCH
    Gamma['axino-N-Z'] = axNZ
    Gamma['axino-N-gamma'] = axNG
    Gamma['axino-C-W'] = axCW
    Gamma['axino-SF-fermion'] = axFSF    
    width = sum(Gamma.values())
    decays = DecayList()
    decays.width = width
    decays.Xfraction = 1.
    if width > 0.:
        decays.addDecay(Decay('axino',['neutralino','radiation'],(axGL+axNH+axCH+axNZ+axNG+axCW+axFSF)/width))
    
 
    return decays