
C....Runs isajet and isared and prints out a minimum of data points to describe the
C....neutralino sigma.v function for the input model.
C....Input (read from terminal):  filename containing the isajet model parameters
C....                              MZ (neutralino1 mass)
C....Output: list of T,sigma.v(T)

      PROGRAM GETISADATA

      IMPLICIT NONE


      DOUBLE PRECISION MZ,T
      CHARACTER*1000 isajet_inputfile,TAG,FILENAME,SIGFILE
      
      LOGICAL DONE,NEWPT
      DOUBLE PRECISION TV(100),SIGV(100),TVNEW(100),SIGVNEW(100)
      DOUBLE PRECISION TNEWV(100)
      DOUBLE PRECISION SIG,TNEW,SIGINT,LININTERP
      INTEGER I,J,NPTS,NPTSNEW,NNEW

      COMMON /SUGMG/ MSS(32),GSS(31),MGUTSS,GGUTSS,AGUTSS,FTGUT,
     $FBGUT,FTAGUT,FNGUT
      REAL MSS,GSS,MGUTSS,GGUTSS,AGUTSS,FTGUT,FBGUT,FTAGUT,FNGUT

      INTEGER IFIRST2
      DOUBLE PRECISION SIG0,SIG1
      COMMON/SFIRSTCALL/SIG0,SIG1,IFIRST2
      SAVE/SFIRSTCALL/

      IFIRST2 = 1
      DO I=1,100
        TV(I) = 0.
        SIGV(I) = 0.
      ENDDO
      NPTSNEW = 0

      READ(*,*) isajet_inputfile
C...Load susy spectrum and compute RGEs (save to file if required in isajet_inputfile)
      OPEN(UNIT=13,FILE=isajet_inputfile,STATUS='OLD')
      CALL PLOTAUX(13) 



C...Generate sigmaV data for neutralino and save to file
      REWIND(13)
      DO I=1,3
        READ(13,*) TAG,FILENAME
        IF (TAG.EQ.'sigmavfile') SIGFILE = FILENAME
      ENDDO  
      CLOSE(13)
      IF (SIGFILE.EQ.'None') RETURN
      OPEN(UNIT=14,FILE=SIGFILE,STATUS='REPLACE')
      MZ = DABS(DBLE(MSS(23)))
C..First sample 11 points, including the extremes:
      DO I=1,11
        TV(I) = (4d-5 + 10.**(DBLE(I)/2.-6.))*MZ
        CALL SIGMAZ(TV(I),MZ,SIGV(I))
        NPTSNEW = NPTSNEW + 1
      ENDDO 
      NPTS = NPTSNEW     
c...Now try to interpolate between the points. Whenever the interpolation fails,
C...include the intermediate points.
      DO I = 1,100
        TNEWV(I) = 0.
      ENDDO
      NNEW = 0
      DONE = .FALSE.
      DO WHILE (.NOT.DONE.AND.NPTS.LT.100)
        NPTSNEW = 0
        DO I=1,NPTS
          SIGVNEW(I) = 0.
          TVNEW(I) = 0.
        ENDDO
        DO I=1,NPTS-1
          NPTSNEW = NPTSNEW + 1
          TVNEW(NPTSNEW) = TV(I)
          SIGVNEW(NPTSNEW) = SIGV(I)          
          TNEW = (TV(I) + TV(I+1))/2.
          NEWPT = .TRUE.
          DO J=1,NNEW
            IF (TNEW.EQ.TNEWV(J)) NEWPT = .FALSE.
          ENDDO
          IF (NEWPT) THEN       !Stores new points and only calls SIGMAZ again if this point is new
            NNEW = NNEW + 1
            TNEWV(NNEW) = TNEW
            CALL SIGMAZ(TNEW,MZ,SIG)
            SIGINT = LININTERP(TV(I),TV(I+1),SIGV(I),SIGV(I+1),TNEW)
            IF (DABS(SIG-SIGINT)/SIG.GT.0.1) THEN   !If interp differs from real value, add intermediate point
              NPTSNEW = NPTSNEW + 1
              SIGVNEW(NPTSNEW) = SIG
              TVNEW(NPTSNEW) = TNEW
            ENDIF
          ENDIF
        ENDDO
        NPTSNEW = NPTSNEW + 1
        TVNEW(NPTSNEW) = TV(NPTS)
        SIGVNEW(NPTSNEW) = SIGV(NPTS)
        IF (NPTSNEW.EQ.NPTS) DONE = .TRUE.   !Stop the loop when the interpolation for all pts are good
        NPTS = NPTSNEW
        DO I=1,NPTS
          TV(I) = TVNEW(I)
          SIGV(I) = SIGVNEW(I)
        ENDDO
      ENDDO

      DO I=1,NPTS
        IF (TV(I).GT.0.) WRITE(14,*) TV(I),SIGV(I)
      ENDDO
      CLOSE(14)

      END PROGRAM


C...Simple linear interpolation function
      DOUBLE PRECISION FUNCTION LININTERP(X1,X2,Y1,Y2,X)

      IMPLICIT NONE

      DOUBLE PRECISION X1,X2,Y1,Y2,X
      DOUBLE PRECISION A,B

      A = (Y2 -Y1)/(X2 - X1)
      B = (X2*Y1 - X1*Y2)/(X2 - X1)

      LININTERP = A*X + B

      END FUNCTION



C...Subroutine to obtain <sigma.v>(T) used to compute the neutralino relic abundance

C...Input
C            T = temperature
C            MZ = neutralino mass
C            SUSY spectrum (parameters set in omegasusy.in)
C            IFIRST2 (if != 166 -> call ISASUSY and ISARED, if == 166 -> assume isasusy and isared were already initialized)

C...Output
C            SIG = <sigma.v>(T) in GeV^-2

C...OBS: For integral convergence issues, SIG is set to SIG(MZ/5) if T > MZ/5
C...FIXES:
C            05/19/2011: Included cut at low temperatures: SIG is set to SIG(MZ/20000) if T < MZ/20000
C            05/19/2011: Now spectrum MUST be initialized in the main program!



      SUBROUTINE SIGMAZ(T,MZ,SIG)

      IMPLICIT NONE

      COMMON /SUGPAS/ XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,NOGOOD,IAL3UN,ITACHY,MHPNEG,IGUTST,ASM3,
     $VUMT,VDMT,ASMTP,ASMSS,M3Q,MHDSQ,MHUSQ,MHDSMG,MHUSMG,MUMG,BMG,
     $FT2Z1,FB2Z1,FL2Z1
      REAL XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,ASM3,VUMT,VDMT,ASMTP,ASMSS,M3Q,MHDSQ,MHUSQ,
     $MHDSMG,MHUSMG,MUMG,BMG,FT2Z1,FB2Z1,FL2Z1
      INTEGER NOGOOD,IAL3UN,ITACHY,MHPNEG,IGUTST
      SAVE /SUGPAS/


      DOUBLE PRECISION AS_MAX,XFINI,XF,COS_MIN,COS_MAX,CS_MIN,SUPPEXP
     &,SUPPEXPMAX
      INTEGER NDIMUSER,NDIMUSER_EFF,ISUM,NST_MAX,
     &NPROC_MIN,NPROC_MAX,NPROC_STEP
      COMMON /CTRL/
     &AS_MAX,XFINI,XF,COS_MIN,COS_MAX,CS_MIN,SUPPEXP, SUPPEXPMAX,
     &NDIMUSER,NDIMUSER_EFF,ISUM,NST_MAX,
     &NPROC_MIN,NPROC_MAX,NPROC_STEP
      SAVE/CTRL/


      INTEGER IFIRST2
      DOUBLE PRECISION SIG0,SIG1
      COMMON/SFIRSTCALL/SIG0,SIG1,IFIRST2
      SAVE/SFIRSTCALL/


      DOUBLE PRECISION T,MZ,SIG,TCMB,RC,RMP,M,TFR,OMH2,PI,GSTAR
      DOUBLE PRECISION FUNC_INT
      EXTERNAL FUNC_INT
      EXTERNAL GSTAR


C...Set 2D integral:
      NDIMUSER=2
      NDIMUSER_EFF=2



C...First compute spectrum and OMGZ if not already initialized
      IF(IFIRST2.NE.166) THEN
        IFIRST2=166
        XF=0.2d0
        SIG0 = FUNC_INT(0)*(2.5682d-9)                                    ! Store limit sigmav value
        IF(NOGOOD.NE.0) THEN
            WRITE(*,*) 'BAD SUSY POINT: STOP'
            GOTO 666
        ENDIF
        XF=5d-5
        SIG1 = FUNC_INT(0)*(2.5682d-9)                                    ! Store limit sigmav value
      ENDIF

C...Set temperature:
      XF=T/MZ
      IF(XF.GT.0.2d0) THEN            ! Freeze sigmav if T >MZ1/5 (to avoid bad integral convergence)
        SIG=SIG0            
      ELSEIF(XF.LT.5d-5) THEN             ! Freeze sigmav if T < MZ1/20000 (to avoid bad integral convergence)      
        SIG=SIG1
      ELSE
C...Compute integral (in pb):
        SIG=FUNC_INT(0)
C...Convert to 1/GeV^2:
        SIG = SIG*(2.5682d-9)
      ENDIF


666      END SUBROUTINE
      SUBROUTINE PLOTAUX(IFL)
C
C     Main subroutine to calculate MSSM input parameters for ISAJET
C     from renormalization group equations and supergravity.
C     All external names are of the form SUxxxx.
C     Must link with block data ALDATA.
C
C     Includes optional link to ISATOOLS, which requires libisared.a.
C     Make this from isared.tar; see the Makefile for instructions.
C
C...##### Reads input from file in unit IFL (plotsugra.in) and writes no output ####
C...!!! Input file must be already opened !!!

C...   10/12/2012: Updated to Isajet v7.83
c...   10/12/2012: When called also saves the running info in TRHPARS
c....04/21/2013: Switched TRHPARS to double precision
c....04/21/2013: Added auxiliary functions PSLAMB and COT


      IMPLICIT NONE
      INTEGER IFL
      COMMON/SSLUN/LOUT,LHEOUT
      INTEGER LOUT,LHEOUT
      SAVE /SSLUN/
C          SUSY parameters
C          AMGLSS               = gluino mass
C          AMULSS               = up-left squark mass
C          AMELSS               = left-selectron mass
C          AMERSS               = right-slepton mass
C          AMNiSS               = sneutrino mass for generation i
C          TWOM1                = Higgsino mass = - mu
C          RV2V1                = ratio v2/v1 of vev's
C          AMTLSS,AMTRSS        = left,right stop masses
C          AMT1SS,AMT2SS        = light,heavy stop masses
C          AMBLSS,AMBRSS        = left,right sbottom masses
C          AMB1SS,AMB2SS        = light,heavy sbottom masses
C          AMLLSS,AMLRSS        = left,right stau masses
C          AML1SS,AML2SS        = light,heavy stau masses
C          AMZiSS               = signed mass of Zi
C          ZMIXSS               = Zi mixing matrix
C          AMWiSS               = signed Wi mass
C          GAMMAL,GAMMAR        = Wi left, right mixing angles
C          AMHL,AMHH,AMHA       = neutral Higgs h0, H0, A0 masses
C          AMHC                 = charged Higgs H+ mass
C          ALFAH                = Higgs mixing angle
C          AAT                  = stop trilinear term
C          THETAT               = stop mixing angle
C          AAB                  = sbottom trilinear term
C          THETAB               = sbottom mixing angle
C          AAL                  = stau trilinear term
C          THETAL               = stau mixing angle
C          AMGVSS               = gravitino mass
C          MTQ                  = top mass at MSUSY
C          MBQ                  = bottom mass at MSUSY
C          MLQ                  = tau mass at MSUSY
C          FBMA                 = b-Yukawa at mA scale
C          VUQ                  = Hu vev at MSUSY
C          VDQ                  = Hd vev at MSUSY
C          SGNM3                = sign of gaugino mass M3
      COMMON/SSPAR/GORGE,AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS(4,4)
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,
     $VUQ,VDQ,SGNM3
      REAL AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,VUQ,VDQ,SGNM3
      LOGICAL GORGE
      REAL AMZISS(4)
      EQUIVALENCE (AMZISS(1),AMZ1SS)
      SAVE /SSPAR/
C          SM ident code definitions. These are standard ISAJET but
C          can be changed.
      INTEGER IDUP,IDDN,IDST,IDCH,IDBT,IDTP
      INTEGER IDNE,IDE,IDNM,IDMU,IDNT,IDTAU
      INTEGER IDGL,IDGM,IDW,IDZ,IDH
      PARAMETER (IDUP=1,IDDN=2,IDST=3,IDCH=4,IDBT=5,IDTP=6)
      PARAMETER (IDNE=11,IDE=12,IDNM=13,IDMU=14,IDNT=15,IDTAU=16)
      PARAMETER (IDGL=9,IDGM=10,IDW=80,IDZ=90,IDH=81)
C          SUSY ident code definitions. They are chosen to be similar
C          to those in versions < 6.50 but may be changed.
      INTEGER ISUPL,ISDNL,ISSTL,ISCHL,ISBT1,ISTP1
      INTEGER ISNEL,ISEL,ISNML,ISMUL,ISNTL,ISTAU1
      INTEGER ISUPR,ISDNR,ISSTR,ISCHR,ISBT2,ISTP2
      INTEGER ISNER,ISER,ISNMR,ISMUR,ISNTR,ISTAU2
      INTEGER ISZ1,ISZ2,ISZ3,ISZ4,ISW1,ISW2,ISGL
      INTEGER ISHL,ISHH,ISHA,ISHC
      INTEGER ISGRAV
      INTEGER IDTAUL,IDTAUR
      PARAMETER (ISUPL=21,ISDNL=22,ISSTL=23,ISCHL=24,ISBT1=25
     &,ISTP1=26)
      PARAMETER (ISNEL=31,ISEL=32,ISNML=33,ISMUL=34,ISNTL=35
     &,ISTAU1=36)
      PARAMETER (ISUPR=41,ISDNR=42,ISSTR=43,ISCHR=44,ISBT2=45
     &,ISTP2=46)
      PARAMETER (ISNER=51,ISER=52,ISNMR=53,ISMUR=54,ISNTR=55
     &,ISTAU2=56)
      PARAMETER (ISGL=29)
      PARAMETER (ISZ1=30,ISZ2=40,ISZ3=50,ISZ4=60,ISW1=39,ISW2=49)
      PARAMETER (ISHL=82,ISHH=83,ISHA=84,ISHC=86)
      PARAMETER (ISGRAV=91)
      PARAMETER (IDTAUL=10016,IDTAUR=20016)
C          Frozen couplings from RG equations:
C     GSS( 1) = g_1        GSS( 2) = g_2        GSS( 3) = g_3
C     GSS( 4) = y_tau      GSS( 5) = y_b        GSS( 6) = y_t
C     GSS( 7) = M_1        GSS( 8) = M_2        GSS( 9) = M_3
C     GSS(10) = A_tau      GSS(11) = A_b        GSS(12) = A_t
C     GSS(13) = M_hd^2     GSS(14) = M_hu^2     GSS(15) = M_er^2
C     GSS(16) = M_el^2     GSS(17) = M_dnr^2    GSS(18) = M_upr^2
C     GSS(19) = M_upl^2    GSS(20) = M_taur^2   GSS(21) = M_taul^2
C     GSS(22) = M_btr^2    GSS(23) = M_tpr^2    GSS(24) = M_tpl^2
C     GSS(25) = mu         GSS(26) = B          GSS(27) = Y_N
C     GSS(28) = M_nr       GSS(29) = A_n        GSS(30) = vdq
C     GSS(31) = vuq
C          Masses:
C     MSS( 1) = glss     MSS( 2) = upl      MSS( 3) = upr
C     MSS( 4) = dnl      MSS( 5) = dnr      MSS( 6) = stl
C     MSS( 7) = str      MSS( 8) = chl      MSS( 9) = chr
C     MSS(10) = b1       MSS(11) = b2       MSS(12) = t1
C     MSS(13) = t2       MSS(14) = nuel     MSS(15) = numl
C     MSS(16) = nutl     MSS(17) = el-      MSS(18) = er-
C     MSS(19) = mul-     MSS(20) = mur-     MSS(21) = tau1
C     MSS(22) = tau2     MSS(23) = z1ss     MSS(24) = z2ss
C     MSS(25) = z3ss     MSS(26) = z4ss     MSS(27) = w1ss
C     MSS(28) = w2ss     MSS(29) = hl0      MSS(30) = hh0
C     MSS(31) = ha0      MSS(32) = h+
C          Unification:
C     MGUTSS  = M_GUT    GGUTSS  = g_GUT    AGUTSS  = alpha_GUT
      COMMON /SUGMG/ MSS(32),GSS(31),MGUTSS,GGUTSS,AGUTSS,FTGUT,
     $FBGUT,FTAGUT,FNGUT
      REAL MSS,GSS,MGUTSS,GGUTSS,AGUTSS,FTGUT,FBGUT,FTAGUT,FNGUT
      SAVE /SUGMG/
C     XSUGIN contains the inputs to SUGRA:
C     XSUGIN(1) = M_0        XSUGIN(2) = M_(1/2)  XSUGIN(3) = A_0
C     XSUGIN(4) = tan(beta)  XSUGIN(5) = sgn(mu)  XSUGIN(6) = M_t
C     XSUGIN(7) = SUG BC scale
C     XGMIN(1) = LAM         XGMIN(2)  = M_MES    XGMIN(3)  = XN5
C     XGMIN(4) = tan(beta)   XGMIN(5)  = sgn(mu)  XGMIN(6) = M_t
C     XGMIN(7) = CGRAV       XGMIN(8)  =RSL       XGMIN(9)  = DEL_HD
C     XGMIN(10)  = DEL_HU    XGMIN(11) = DY       XGMIN(12) = N5_1
C     XGMIN(13)  = N5_2      XGMIN(14) = N5_3
C     XNRIN(1) = M_N3        XNRIN(2) = M_MAJ     XNRIN(3) = ANSS 
C     XNRIN(4) = M_N3SS
C     XISAIN contains the MSSMi inputs in natural order.
      COMMON /SUGXIN/ XISAIN(24),XSUGIN(7),XGMIN(14),XNRIN(4),
     $XAMIN(11)
      REAL XISAIN,XSUGIN,XGMIN,XNRIN,XAMIN
      SAVE /SUGXIN/
      COMMON /SUGPAS/ XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,NOGOOD,IAL3UN,ITACHY,MHPNEG,MHLNEG,MHCNEG,
     $MSQNEG,IGUTST,ASM3,
     $VUMT,VDMT,ASMTP,ASMSS,M3Q,MHDSQ,MHUSQ,MHDSMG,MHUSMG,MUMG,BMG,
     $FT2Z1,FB2Z1,FL2Z1,SIGDMX,SIGUMX,C5MAX,C5MAXV(20)
      REAL XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,ASM3,VUMT,VDMT,ASMTP,ASMSS,M3Q,MHDSQ,MHUSQ,
     $MHDSMG,MHUSMG,MUMG,BMG,FT2Z1,FB2Z1,FL2Z1,
     $SIGDMX,SIGUMX,C5MAX,C5MAXV
      INTEGER NOGOOD,IAL3UN,ITACHY,MHPNEG,MHLNEG,MHCNEG,
     $MSQNEG,IGUTST
      SAVE /SUGPAS/
C     XNUSUG contains non-universal GUT scale soft terms for SUGRA:
C     XNUSUG(1)=M1 XNUSUG(2)=M2 XNUSUG(3)=M3
C     XNUSUG(4)=A_tau XNUSUG(5)=A_b XNUSUG(6)=A_t
C     XNUSUG(7)=m_Hd XNUSUG(8)=m_Hu XNUSUG(9)=m_eR XNUSUG(10)=m_eL
C     XNUSUG(11)=m_dR XNUSUG(12)=m_uR XNUSUG(13)=m_uL XNUSUG(14)=m_lR
C     XNUSUG(15)=m_lL XNUSUG(16)=m_bR XNUSUG(17)=m_tR XNUSUG(18)=m_tL
C     XNUSUG(19)=mu(Q) XNUSUG(20)=mA(Q)
      COMMON /SUGNU/ XNUSUG(20),INUHM
      REAL XNUSUG
      INTEGER INUHM
      SAVE /SUGNU/
C          ISAPW1 is used to check whether ALDATA is loaded
      COMMON/ISAPW/ISAPW1
      CHARACTER*30 ISAPW1
      SAVE /ISAPW/
      CHARACTER*80 FNAME,FNLHA,FNWIG
      LOGICAL GOLHA,GOWIG
      INTEGER ILHA,IWIG,IMHL,IMHC,IMSQ
      REAL M0,MHF,A0,TANB,SGNMU,MT,XLAMGM,XMESGM,XN5GM,AMPL,XCMGV
      INTEGER NSTEP,IMODEL,INUSUG,IMODIN,ISATLS,IRGEFL
      INTEGER K,NOUT,IALLOW,IITEST,J,II
      CHARACTER*40 VERSN,VISAJE
      PARAMETER (NOUT=33)
      INTEGER IDOUT(NOUT)
      CHARACTER*30 ISAPW2
      SAVE ISAPW2
C
C          Isatools common blocks and variables
C
      COMMON/SUGRED/OMGH2,SIGMA,XFREEZ,NSTEPS,FFF_V
      REAL OMGH2,SIGMA,XFREEZ,FFF_V
      INTEGER NSTEPS
      REAL ALEMIGM2,BFBSG,ALEMI
      COMMON/SUGRES/SIGMA0PROT,SIGMA0NEUT,SIGMASPROT,SIGMASNEUT
      REAL*8 SIGMA0PROT,SIGMA0NEUT,SIGMASPROT,SIGMASNEUT
      SAVE/SUGRES/
C-FP  INTEGER INUHM
      COMMON/RGEFNM/FNRGE
      CHARACTER*128 FNRGE
      REAL*8 DAMU,DBFBSG
      REAL BRBS,BRBD
      REAL R,BFBTN
      INTEGER IRED,IRES,IAMU,IBSG,IBLL,IBTN
C
      DATA IDOUT/
     $IDTP,ISGL,ISUPL,ISDNL,ISSTL,ISCHL,ISBT1,ISTP1,ISUPR,ISDNR,
     $ISSTR,ISCHR,ISBT2,ISTP2,ISEL,ISMUL,ISTAU1,ISNEL,ISNML,ISNTL,
     $ISER,ISMUR,ISTAU2,ISZ1,ISZ2,ISZ3,ISZ4,ISW1,ISW2,
     $ISHL,ISHH,ISHA,ISHC/
      DATA AMPL/2.4E18/
C          ISAPW2 is used to check whether ALDATA is loaded
      DATA ISAPW2/'ALDATA REQUIRED BY FORTRAN G,H'/
      DATA ILHA/11/,IWIG/12/

      DOUBLE PRECISION TRHM1V,TRHM2V,TRHM3V,TRHG1V,TRHG2V,TRHG3V,TQV
      COMMON/TRHPARS/ TRHM1V(100000),TRHM2V(100000),TRHM3V(100000),
     &TRHG1V(100000),TRHG2V(100000),TRHG3V(100000),TQV(100000)
      CHARACTER*1000 SLHAFILE,RGEFILE,TAG,FILENAME
      INTEGER IRGE,I
      LOGICAL PRINTRGE


      INTEGER IDLSP

C
C          Initialize
C


      REWIND(IFL)

      IF(ISAPW1.NE.ISAPW2) THEN
        PRINT*, ' ERROR: BLOCK DATA ALDATA HAS NOT BEEN LOADED.'
        PRINT*, ' ISAJET CANNOT RUN WITHOUT IT.'
        PRINT*, ' PLEASE READ THE FINE MANUAL FOR ISAJET.'
        STOP99
      ENDIF
C
      LOUT=1
      NSTEP=1000
      XNRIN(2)=1.E20
      XSUGIN(7)=0.
      INUHM=0
C
C          Open files
C
      DO I=1,3
        READ(IFL,*) TAG,FILENAME
        IF (TAG.EQ.'rgefile') THEN
          RGEFILE = FILENAME
        ELSE IF (TAG.EQ.'slhafile') THEN
          SLHAFILE = FILENAME
        ELSE IF (TAG.EQ.'sigmavfile') THEN
          IF (FILENAME.EQ.'None') ISATLS = -1   !Do not run isatools if sigmav data is not needed
        ENDIF
      ENDDO
      IF (RGEFILE.EQ.'None') THEN
        PRINTRGE=.FALSE.      
      ELSE
        PRINTRGE=.TRUE.      ! Prints RGE evolution of gauge couplings and gaugino masses to file
        OPEN(IRGE,FILE=RGEFILE,STATUS='REPLACE',FORM='FORMATTED')
      ENDIF
      IF (SLHAFILE.EQ.'None') THEN
        GOLHA=.FALSE.      ! No les houche file
      ELSE
        GOLHA=.TRUE. 
        OPEN(ILHA,FILE=SLHAFILE,STATUS='REPLACE',FORM='FORMATTED')
      ENDIF
      GOWIG=.FALSE.      ! No herwig file
      READ(IFL,*)  IMODIN
      IMODEL=IMODIN
      IF (IMODEL.EQ.4) THEN
        IAL3UN=1
        IMODEL=1
      END IF
      IF (IMODEL.EQ.1.OR.IMODEL.EQ.3.OR.IMODEL.EQ.6) THEN
        READ(IFL,*)  M0,MHF,A0,TANB,SGNMU,MT
        IF (IMODEL.EQ.6) THEN
          IMODEL=1
          READ(IFL,*)  XNRIN(1),XNRIN(2),XNRIN(3),XNRIN(4)
          GO TO 15
        END IF
        IF (IMODEL.EQ.3) THEN
10        IMODEL=1
          READ(IFL,*)  INUSUG
          IF (INUSUG.EQ.0) THEN
            GO TO 15
          ELSE IF (INUSUG.EQ.1) THEN
C            PRINT*,'Enter GUT scale M_1, M_2, M_3:'
            READ(IFL,*)  XNUSUG(1),XNUSUG(2),XNUSUG(3)
C            IF (XNUSUG(3).LE.0.) THEN
C            PRINT*, ' NEGATIVE M_3 IS NOT ALLOWED'
C            STOP 99
C            END IF
          ELSE IF (INUSUG.EQ.2) THEN
C            PRINT*,'Enter GUT scale A_t, A_b, A_tau:'
            READ(IFL,*)  XNUSUG(6),XNUSUG(5),XNUSUG(4)
          ELSE IF (INUSUG.EQ.3) THEN
C            PRINT*,'Enter GUT scale m_Hd, m_Hu:'
            READ(IFL,*)  XNUSUG(7),XNUSUG(8)
          ELSE IF (INUSUG.EQ.4) THEN
C            PRINT*,'Enter GUT scale M(ul), M(dr), M(ur), M(el), M(er):'
            READ(IFL,*)  XNUSUG(13),XNUSUG(11),XNUSUG(12),XNUSUG(10),
     $      XNUSUG(9)
          ELSE IF (INUSUG.EQ.5) THEN
C            PRINT*,'Enter GUT scale M(tl), M(br), M(tr), M(Ll), M(Lr):'
            READ(IFL,*)  XNUSUG(18),XNUSUG(16),XNUSUG(17),XNUSUG(15),
     $      XNUSUG(14)
          ELSE IF (INUSUG.EQ.6) THEN
C            PRINT*,' ENTER M(nu_3), M_Majorana, A_N, M(NRSS)'
            READ(IFL,*)  XNRIN(1),XNRIN(2),XNRIN(3),XNRIN(4)
          ELSE IF (INUSUG.EQ.7) THEN
C            PRINT*,' ENTER Q_max high scale for SUSY BCs'
            READ(IFL,*)  XSUGIN(7)
          ELSE IF (INUSUG.EQ.8) THEN
C            PRINT*,' ENTER mu(Q), mA(Q)'
            READ(IFL,*)  XNUSUG(19),XNUSUG(20)
            MU=XNUSUG(19)
            AMHA=XNUSUG(20)
            TWOM1=-MU
            INUHM=1
          END IF
          GO TO 10
        END IF
      ELSE IF (IMODEL.EQ.2.OR.IMODEL.EQ.5) THEN
C          PRINT*,'ENTER Lambda, M_mes, N_5, tan(beta), sgn(mu), ',
C     $    'M_t, C_gv:'
          READ(IFL,*)  M0,MHF,A0,TANB,SGNMU,MT,XCMGV
          XGMIN(7)=XCMGV
          XGMIN(8)=1.
          AMGVSS=M0*MHF*XCMGV/SQRT(3.)/AMPL
          IF (IMODEL.EQ.5) THEN
            IMODEL=2
C            PRINT*,'Rsl = factor multiplying gaugino masses at M_mes'
C            PRINT*,'dmH_d^2, dmH_u^2 = Higgs mass**2 shifts at M_mes'
C            PRINT*,'d_Y = mass**2 shifts proportional to Y at M_mes'
C            PRINT*,'n5_1,n5_2,n5_3 = n5 values for U(1),SU(2),SU(3)'
C            PRINT*,'ENTER Rsl, dmH_d^2, dmH_u^2, d_Y, n5_1, n5_2, n5_3'
            READ(IFL,*)  XGMIN(8),XGMIN(9),XGMIN(10),XGMIN(11),
     $      XGMIN(12),
     $      XGMIN(13),XGMIN(14)
            END IF
      ELSE IF (IMODEL.EQ.7) THEN
C        PRINT*,'ENTER M_0, M_(3/2), tan(beta), sgn(mu), M_t:'
        READ(IFL,*)  M0,MHF,TANB,SGNMU,MT
        A0=0.
        DO 101 II=1,7
101     XAMIN(II)=1.
      ELSE IF (IMODEL.EQ.8) THEN
C        PRINT*,'ENTER M_0, M_(3/2), tan(beta), sgn(mu), M_t:'
        READ(IFL,*)  M0,MHF,TANB,SGNMU,MT
        A0=0.
C        PRINT*,'ENTER cQ, cD, cU, cL, cE, cHd, cHu:'
        READ(IFL,*)  
     $XAMIN(1),XAMIN(2),XAMIN(3),XAMIN(4),XAMIN(5),XAMIN(6),XAMIN(7)
        IMODEL=7
      ELSE IF (IMODEL.EQ.9) THEN
C        PRINT*,'ENTER alpha, M_(3/2), tan(beta), sgn(mu), M_t:'
        READ(IFL,*)  M0,MHF,TANB,SGNMU,MT
        A0=0.
C          Set defaults
        DO 102 II=1,7
          XAMIN(II)=0
102     CONTINUE
        DO 103 II=8,10
          XAMIN(II)=1
103     CONTINUE
C        PRINT*,'ENTER moduli weights nQ, nD, nU, nL, nE, nHd, nHu ',
c     $  '[/ for all 0]:'
        READ(IFL,*)  
     $XAMIN(1),XAMIN(2),XAMIN(3),XAMIN(4),XAMIN(5),XAMIN(6),XAMIN(7)
C        PRINT*,'ENTER moduli parameters L1, L2, L3 [/ for all 1]:'
        READ(IFL,*)  XAMIN(8),XAMIN(9),XAMIN(10)
      ELSE IF (IMODEL.EQ.10) THEN
C        PRINT*,'ENTER alpha, M_(3/2), tan(beta), sgn(mu), M_t:'
        READ(IFL,*)  XAMIN(11),MHF,TANB,SGNMU,MT
        M0=0.
        A0=0.
      ELSE
        PRINT*,'Invalid model choice.'
        STOP99
      END IF

C
C          Solve RG equations
C
15    CALL MODSUGRA(M0,MHF,A0,TANB,SGNMU,MT,IMODEL)
c15    CALL SUGRA(M0,MHF,A0,TANB,SGNMU,MT,IMODEL)
C
C          Print results
C
      
      VERSN=VISAJE()
c      WRITE(LOUT,20) VERSN
20    FORMAT(' ',44('*')/' *',42X,'*'/
     $  ' * ',A40,' *'/
     $  ' *',42X,'*'/' ',44('*')/)
      IF (NOGOOD.EQ.1) THEN
        PRINT*, 'BAD POINT: TACHYONIC PARTICLES!'
c        WRITE(LOUT,*) 'BAD POINT: TACHYONIC PARTICLES!'
      ELSE IF (NOGOOD.EQ.2) THEN
        PRINT*, 'BAD POINT: NO EW SYMMETRY BREAKING!'
c        WRITE(LOUT,*) 'BAD POINT: NO EW SYMMETRY BREAKING!'
      ELSE IF (NOGOOD.EQ.3) THEN
        PRINT*, 'BAD POINT: M(H_P)^2<0!'
c        WRITE(LOUT,*) 'BAD POINT: M(H_P)^2<0!'
      ELSE IF (NOGOOD.EQ.4) THEN
        PRINT*, 'BAD POINT: YUKAWA>10!'
c        WRITE(LOUT,*) 'BAD POINT: YUKAWA>10!'
      ELSE IF (NOGOOD.EQ.5.AND.IMODEL.EQ.1) THEN
        PRINT*, 'SUGRA BAD POINT: Z1SS NOT LSP!'
c        WRITE(LOUT,*) 'SUGRA BAD POINT: Z1SS NOT LSP!'
      ELSE IF (NOGOOD.EQ.7) THEN
        PRINT*, 'BAD POINT: XT EWSB BAD!'
c        WRITE(LOUT,*) 'BAD POINT: XT EWSB BAD!'
      ELSE IF (NOGOOD.EQ.8) THEN
        PRINT*, 'BAD POINT: MHL^2<0!'
c        WRITE(LOUT,*) 'BAD POINT: MHL^2<0!'
      ELSE IF (NOGOOD.EQ.-1) THEN
        PRINT*, 'BAD POINT: NO RGE SOLUTION FOUND'
c        WRITE(LOUT,*) 'BAD POINT: NO RGE SOLUTION FOUND'
      END IF
      IF (MHPNEG.EQ.1) THEN
        PRINT*, 'BAD POINT: M(H_P)^2<0!!'
c        WRITE(LOUT,*) 'BAD POINT: M(H_P)^2<0!!'
        NOGOOD=3
      END IF



c      IF(NOGOOD.NE.0) STOP99
      IF(ITACHY.NE.0) THEN
c        WRITE(LOUT,*) 'WARNING: TACHYONIC SLEPTONS AT GUT SCALE'
c        WRITE(LOUT,*) '         POINT MAY BE INVALID'
      ENDIF
CCC      IF (IGUTST.EQ.1) THEN
CCC        PRINT*, 'WARNING: GUT INSTABILITY IN NUHM MODEL'
CCC      END IF
C
C          Print selected model and results
C
c      IF(IMODIN.EQ.1) WRITE(LOUT,1001)
1001  FORMAT(//' Minimal supergravity (mSUGRA) model:'/)
c      IF(IMODIN.EQ.2) WRITE(LOUT,1002)
1002  FORMAT(//' Minimal gauge mediated (GMSB) model:'/)
c      IF(IMODIN.EQ.3) WRITE(LOUT,1003)
1003  FORMAT(//' Non-universal supergravity model:'/)
c      IF(IMODIN.EQ.4) WRITE(LOUT,1004)
1004  FORMAT(//' Supergravity model with truly unified couplings:'/)
c      IF(IMODIN.EQ.5) WRITE(LOUT,1005)
1005  FORMAT(//' Non-minimal gauge mediated (GMSB) model:'/)
c      IF(IMODIN.EQ.6) WRITE(LOUT,1006)
1006  FORMAT(//' Supergravity model with right-handed neutrinos:'/)
c      IF(IMODIN.EQ.7) WRITE(LOUT,1007)
1007  FORMAT(//' Anomaly-mediated SUSY breaking model:'/)
c      IF(IMODIN.EQ.8) WRITE(LOUT,1008)
1008  FORMAT(//' Non-minimal anomaly-mediated SUSY breaking model:'/)
c      IF(IMODIN.EQ.9) WRITE(LOUT,1009)
1009  FORMAT(//' Mixed modulus-AMSB SUSY breaking model:'/)
c      IF(IMODIN.EQ.10) WRITE(LOUT,1010)
1010  FORMAT(//' Hypercharged-AMSB SUSY breaking model:'/)
C
C          Calculate all masses and decay modes
C
 
C
C          Calculate all masses and decay modes
C
        CALL SSMSSM(XISAIN(1),XISAIN(2),XISAIN(3),
     $ XISAIN(4),XISAIN(5),XISAIN(6),XISAIN(7),XISAIN(8),XISAIN(9),
     $ XISAIN(10),XISAIN(11),XISAIN(12),XISAIN(13),XISAIN(14),
     $ XISAIN(15),XISAIN(16),XISAIN(17),XISAIN(18),XISAIN(19),
     $ XISAIN(20),XISAIN(21),XISAIN(22),XISAIN(23),XISAIN(24),
     $ MT,IALLOW,IMODEL,IMHL,IMHC,IMSQ)


      CALL FINDLSP(IDLSP)
      IF(IDLSP.NE.23) THEN
        NOGOOD=5
        PRINT*, 'BAD POINT: Z1SS NOT LSP!',MHF,IDLSP,MSS(23),
     &MSS(IDLSP)
      ENDIF

C
C          Execute Isatools
C

      IRED=0
      IRES=0
      IAMU=0
      IBSG=0
      IBLL=0
      IBTN=0
C      PRINT*,'Run Isatools? Choose 2=all, 1=some, 0=none:'
      IF(ISATLS.EQ.-1) THEN      ! Ignore isatools option
         ISATLS=0
         READ(IFL,*)
      ELSE
         READ(IFL,*)  ISATLS
      ENDIF
      IF(ISATLS.EQ.2) THEN
        IRED=1
        IRES=1
        IAMU=1
        IBSG=1
        IBLL=1
      ELSE IF (ISATLS.EQ.1) THEN
C        PRINT*,'Select desired ISATools packages'
C        PRINT*,'Neutralino Relic Density [1/0]:'
        READ(IFL,*)  IRED
C        PRINT*,'Neutralino DD rates [1/0]:'
        READ(IFL,*)  IRES
C        PRINT*,'Muon (g-2)/2 [1/0]:'
        READ(IFL,*)  IAMU
C        PRINT*,'b->s gamma branching fraction [1/0]:'
        READ(IFL,*)  IBSG
C        PRINT*,'B_s->ll branching fractions [1/0]:'
        READ(IFL,*)  IBLL
      ENDIF


      DAMU=0.
      BFBSG=0.
      OMGH2=0.
      BRBS=0.
      BRBD=0.
      BFBTN=0.
      SIGMA0PROT=0.
      SIGMA0NEUT=0.
      SIGMASPROT=0.
      SIGMASNEUT=0.


C...Don't run IsaTools if point is not good:
      IF(NOGOOD.NE.0) GOTO 213


C          g_mu - 2
      IF (IAMU.EQ.1) THEN
        ALEMI=128.
      CALL ISAAMU(RV2V1,ALEMI,GAMMAL,GAMMAR,TWOM1,AAL,SQRT(GSS(16)),
     $SQRT(GSS(15)),AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,AMW1SS,AMW2SS,AMN2SS,
     $              ZMIXSS,0,DAMU)
      ENDIF
C          B -> s gamma
      IF (IBSG.EQ.1) THEN
        CALL ISABSG(IMODEL,M0,MHF,A0,DBFBSG,0)
        BFBSG=SNGL(DBFBSG)
      ENDIF
C          Bu -> tau+nu_tau
      IF (IBTN.EQ.1) THEN
        CALL ISABTN(TANB,MSS(32),MU,MSS(1),MSS(10),MSS(11),R,BFBTN)
      ENDIF
C          Dark matter cross sections
      IF (IRES.EQ.1) THEN
        CALL ISARES(0)
      ENDIF
C          Relic density calculation -- requires isared.tar
      IF (IRED.EQ.1) THEN
        CALL ISARED(0)
      ENDIF
C          B_s -> mu mu (BRBS) and B -> tau tau (BRBD)
      IF (IBLL.EQ.1) THEN
        CALL ISABMM(MT,TANB,MSS(29),MSS(30),MSS(31),MSS(1),
     $              MSS(10),MSS(11),GSS(11),THETAB,
     $              MSS(12),MSS(13),GSS(12),THETAT,
     $              MU,GSS(8),ALFAH,MSS(6),MSS(8),BRBS,BRBD)
      ENDIF
C
c      CALL SUGPRT(IMODEL,IMODIN)
C
C          Test parameters
C
      IF(IALLOW.NE.0) THEN
c        WRITE(LOUT,2001)
2001    FORMAT(//' MSSM WARNING: Z1SS IS NOT LSP')
      ENDIF
C
      CALL SSTEST(IALLOW)
      IITEST=IALLOW/2
      IF(MOD(IITEST,2).NE.0) THEN
c        WRITE(LOUT,2002)
2002    FORMAT(' MSSM WARNING: Z -> Z1SS Z1SS EXCEEDS BOUND')
      ENDIF
      IITEST=IITEST/2
      IF(MOD(IITEST,2).NE.0) THEN
c        WRITE(LOUT,2004)
2004    FORMAT(' MSSM WARNING: Z -> CHARGINOS ALLOWED')
      ENDIF
      IITEST=IITEST/2
      IF(MOD(IITEST,2).NE.0) THEN
c        WRITE(LOUT,2008)
2008    FORMAT(' MSSM WARNING: Z -> Z1SS Z2SS TOO BIG')
      ENDIF
      IITEST=IITEST/2
      IF(MOD(IITEST,2).NE.0) THEN
c        WRITE(LOUT,2016)
2016    FORMAT(' MSSM WARNING: Z -> SQUARKS, SLEPTONS ALLOWED')
      ENDIF
      IITEST=IITEST/2
      IF(MOD(IITEST,2).NE.0) THEN
c        WRITE(LOUT,2032)
2032    FORMAT(' MSSM WARNING: Z -> Z* HL0 EXCEEDS BOUND')
      ENDIF
      IITEST=IITEST/2
      IF(MOD(IITEST,2).NE.0) THEN
c        WRITE(LOUT,2064)
2064    FORMAT(' MSSM WARNING: Z -> HL0 HA0 ALLOWED')
      ENDIF
C
C          Isatools output
C
      BFBSG=1.0E-4*BFBSG
c      WRITE(LOUT,3500) DAMU,BFBSG,OMGH2,FFF_V,BRBS,BRBD,BFBTN,
c     $SIGMA0PROT,SIGMA0NEUT,SIGMASPROT,SIGMASNEUT
3500  FORMAT(//' Output from ISATOOLS:'/
     $' Delta a_mu       = ',E12.5,'    BF(b->s gamma)    = ',E12.5/
     $' Omega h^2        = ',E12.5,'    <s.v>(v=0)[cm^3/s]= ',E12.5/
     $' BF(Bs -> mu mu)  = ',E12.5,'    BF(B -> tau tau)  = ',E12.5/
     $' BF(Bu -> tau nu) = ',E12.5/
     $' LSP-nucleon spin independent (SI) and dependent (SD) sigmas:'/
     $' sigma(p,SI)[pb]  = ',E12.5,'    sigma(n,SI)[pb]   = ',E12.5/
     $' sigma(p,SD)[pb]  = ',E12.5,'    sigma(n,SD)[pb]   = ',E12.5)
c      WRITE(LOUT,3600)
3600  FORMAT(//' ISASUSY decay modes:'/
     $' Parent --> daughters',18X,'Width',10X,'Branching ratio'/)
C          Write all modes
c      DO 200 J=1,NOUT
c        CALL SSPRT(IDOUT(J))
c200   CONTINUE
C
C          Make optional output files
C

      IF(PRINTRGE) THEN
        DO I=1,100000
          TRHG1V(I) = TRHG1V(I)*SQRT(3./5.)
          IF (TQV(I).NE.0.) WRITE(IRGE,'(10(E10.3,1X))') TQV(I)
     &,TRHM1V(I),TRHM2V(I),TRHM3V(I),TRHG1V(I),TRHG2V(I),TRHG3V(I)
        ENDDO
      ENDIF

      IF(GOLHA) THEN
        CALL ISALHA(ILHA,IMODEL,IMODIN,MT)
        DO 210 J=1,NOUT
          CALL ISALHD(ILHA,IDOUT(J),J,NOUT)
210     CONTINUE
      ENDIF
      IF(GOWIG) THEN
        CALL ISAWIG(IWIG,0,MT,GSS(8),SQRT(GSS(16)),SQRT(GSS(15)),
     &              SQRT(GSS(19)),SQRT(GSS(18)),SQRT(GSS(17)),
     &              SQRT(GSS(16)),SQRT(GSS(15)),SQRT(GSS(19)),
     &              SQRT(GSS(18)),SQRT(GSS(17)),SQRT(GSS(21)),
     &              SQRT(GSS(20)),SQRT(GSS(24)),SQRT(GSS(23)),
     &              SQRT(GSS(22)))
      ENDIF
C
      FNRGE=''
      IRGEFL=0
c      PRINT*,'To run RGEFLAV, enter filename Prefix [/ for none]:'
c      PRINT*,'RGEFLAV will open file Prefix.rgein, and print to files'
c      PRINT*,'Prefix.weakout, Prefix.gutout, Prefix.sqm2u, Prefix.sqm2d'
c      READ*,FNRGE
c      IF(FNRGE.NE.'') THEN
c        IRGEFL=1
c      ENDIF
c      IF(IRGEFL.EQ.1) THEN
c        CALL RGEFLAV(M0,MHF,A0,TANB,SGNMU,MT)
c      END IF
C
213   REWIND(IFL)
C 
c      STOP
      END

c....Modified version of SUGRA SUBROUTINE. Calls MODSUGRGE instead of SUGRGE (only change)
C--------------------------------------------------------------------
      SUBROUTINE MODSUGRA(M0,MHF,A0,TANB,SGNMU,MT,IMODEL)
C--------------------------------------------------------------------
C
C     Calculate supergravity spectra for ISAJET using as inputs
C     M0    = M_0       = common scalar mass at GUT scale
C     MHF   = M_(1/2)   = common gaugino mass at GUT scale
C     A0    = A_0       = trilinear soft breaking parameter at GUT scale
C     TANB  = tan(beta) = ratio of vacuum expectation values v_1/v_2
C     SGNMU = sgn(mu)   = +-1 = sign of Higgsino mass term
C     MT    = M_t       = mass of t quark
C     M0    = Lambda    = ratio of vevs <F>/<S>
C     MHF   = M_Mes     = messenger scale
C     A0    = n_5       = number of messenger fields
C     IMODEL            = 1 for SUGRA model
C                       = 2 for GMSB model
C                       = 7 for AMSB model
C
C     Uses Runge-Kutta method to integrate RGE's from M_Z to M_GUT
C     and back, putting in correct thresholds. For the first iteration
C     only the first 6 couplings are included and a common threshold
C     is used.
C
C     See /SUGMG/ for definitions of couplings and masses.
C
C     Ver. 7.64: Use different convergence cuts for Higgs-related
C                soft couplings.
C     Ver. 7.70 Implement double precision RGE running
C
      IMPLICIT NONE
      COMMON/SSLUN/LOUT,LHEOUT
      INTEGER LOUT,LHEOUT
      SAVE /SSLUN/
C          SUSY parameters
C          AMGLSS               = gluino mass
C          AMULSS               = up-left squark mass
C          AMELSS               = left-selectron mass
C          AMERSS               = right-slepton mass
C          AMNiSS               = sneutrino mass for generation i
C          TWOM1                = Higgsino mass = - mu
C          RV2V1                = ratio v2/v1 of vev's
C          AMTLSS,AMTRSS        = left,right stop masses
C          AMT1SS,AMT2SS        = light,heavy stop masses
C          AMBLSS,AMBRSS        = left,right sbottom masses
C          AMB1SS,AMB2SS        = light,heavy sbottom masses
C          AMLLSS,AMLRSS        = left,right stau masses
C          AML1SS,AML2SS        = light,heavy stau masses
C          AMZiSS               = signed mass of Zi
C          ZMIXSS               = Zi mixing matrix
C          AMWiSS               = signed Wi mass
C          GAMMAL,GAMMAR        = Wi left, right mixing angles
C          AMHL,AMHH,AMHA       = neutral Higgs h0, H0, A0 masses
C          AMHC                 = charged Higgs H+ mass
C          ALFAH                = Higgs mixing angle
C          AAT                  = stop trilinear term
C          THETAT               = stop mixing angle
C          AAB                  = sbottom trilinear term
C          THETAB               = sbottom mixing angle
C          AAL                  = stau trilinear term
C          THETAL               = stau mixing angle
C          AMGVSS               = gravitino mass
C          MTQ                  = top mass at MSUSY
C          MBQ                  = bottom mass at MSUSY
C          MLQ                  = tau mass at MSUSY
C          FBMA                 = b-Yukawa at mA scale
C          VUQ                  = Hu vev at MSUSY
C          VDQ                  = Hd vev at MSUSY
C          SGNM3                = sign of gaugino mass M3
      COMMON/SSPAR/GORGE,AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS(4,4)
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,
     $VUQ,VDQ,SGNM3
      REAL AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,VUQ,VDQ,SGNM3
      LOGICAL GORGE
      REAL AMZISS(4)
      EQUIVALENCE (AMZISS(1),AMZ1SS)
      SAVE /SSPAR/
C          Standard model parameters
C          AMUP,...,AMTP        = quark masses
C          AME,AMMU,AMTAU       = lepton masses
C          AMW,AMZ              = W,Z masses
C          GAMW,GAMZ            = W,Z widths
C          ALFAEM,SN2THW,ALFA3  = SM couplings
C          ALQCD4               = 4 flavor lambda
      COMMON/SSSM/AMUP,AMDN,AMST,AMCH,AMBT,AMTP,AME,AMMU,AMTAU
     $,AMW,AMZ,GAMW,GAMZ,ALFAEM,SN2THW,ALFA2,ALFA3,ALQCD4
      REAL AMUP,AMDN,AMST,AMCH,AMBT,AMTP,AME,AMMU,AMTAU
     $,AMW,AMZ,GAMW,GAMZ,ALFAEM,SN2THW,ALFA2,ALFA3,ALQCD4
      SAVE /SSSM/
C     XSUGIN contains the inputs to SUGRA:
C     XSUGIN(1) = M_0        XSUGIN(2) = M_(1/2)  XSUGIN(3) = A_0
C     XSUGIN(4) = tan(beta)  XSUGIN(5) = sgn(mu)  XSUGIN(6) = M_t
C     XSUGIN(7) = SUG BC scale
C     XGMIN(1) = LAM         XGMIN(2)  = M_MES    XGMIN(3)  = XN5
C     XGMIN(4) = tan(beta)   XGMIN(5)  = sgn(mu)  XGMIN(6) = M_t
C     XGMIN(7) = CGRAV       XGMIN(8)  =RSL       XGMIN(9)  = DEL_HD
C     XGMIN(10)  = DEL_HU    XGMIN(11) = DY       XGMIN(12) = N5_1
C     XGMIN(13)  = N5_2      XGMIN(14) = N5_3
C     XNRIN(1) = M_N3        XNRIN(2) = M_MAJ     XNRIN(3) = ANSS
C     XNRIN(4) = M_N3SS
C     XISAIN contains the MSSMi inputs in natural order.
      COMMON /SUGXIN/ XISAIN(24),XSUGIN(7),XGMIN(14),XNRIN(4),
     $XAMIN(11)
      REAL XISAIN,XSUGIN,XGMIN,XNRIN,XAMIN
      SAVE /SUGXIN/
C          Frozen couplings from RG equations:
C     GSS( 1) = g_1        GSS( 2) = g_2        GSS( 3) = g_3
C     GSS( 4) = y_tau      GSS( 5) = y_b        GSS( 6) = y_t
C     GSS( 7) = M_1        GSS( 8) = M_2        GSS( 9) = M_3
C     GSS(10) = A_tau      GSS(11) = A_b        GSS(12) = A_t
C     GSS(13) = M_hd^2     GSS(14) = M_hu^2     GSS(15) = M_er^2
C     GSS(16) = M_el^2     GSS(17) = M_dnr^2    GSS(18) = M_upr^2
C     GSS(19) = M_upl^2    GSS(20) = M_taur^2   GSS(21) = M_taul^2
C     GSS(22) = M_btr^2    GSS(23) = M_tpr^2    GSS(24) = M_tpl^2
C     GSS(25) = mu         GSS(26) = B          GSS(27) = Y_N
C     GSS(28) = M_nr       GSS(29) = A_n        GSS(30) = vdq
C     GSS(31) = vuq
C          Masses:
C     MSS( 1) = glss     MSS( 2) = upl      MSS( 3) = upr
C     MSS( 4) = dnl      MSS( 5) = dnr      MSS( 6) = stl
C     MSS( 7) = str      MSS( 8) = chl      MSS( 9) = chr
C     MSS(10) = b1       MSS(11) = b2       MSS(12) = t1
C     MSS(13) = t2       MSS(14) = nuel     MSS(15) = numl
C     MSS(16) = nutl     MSS(17) = el-      MSS(18) = er-
C     MSS(19) = mul-     MSS(20) = mur-     MSS(21) = tau1
C     MSS(22) = tau2     MSS(23) = z1ss     MSS(24) = z2ss
C     MSS(25) = z3ss     MSS(26) = z4ss     MSS(27) = w1ss
C     MSS(28) = w2ss     MSS(29) = hl0      MSS(30) = hh0
C     MSS(31) = ha0      MSS(32) = h+
C          Unification:
C     MGUTSS  = M_GUT    GGUTSS  = g_GUT    AGUTSS  = alpha_GUT
      COMMON /SUGMG/ MSS(32),GSS(31),MGUTSS,GGUTSS,AGUTSS,FTGUT,
     $FBGUT,FTAGUT,FNGUT
      REAL MSS,GSS,MGUTSS,GGUTSS,AGUTSS,FTGUT,FBGUT,FTAGUT,FNGUT
      SAVE /SUGMG/
      COMMON /SUGPAS/ XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,NOGOOD,IAL3UN,ITACHY,MHPNEG,MHLNEG,MHCNEG,
     $MSQNEG,IGUTST,ASM3,
     $VUMT,VDMT,ASMTP,ASMSS,M3Q,MHDSQ,MHUSQ,MHDSMG,MHUSMG,MUMG,BMG,
     $FT2Z1,FB2Z1,FL2Z1,SIGDMX,SIGUMX,C5MAX,C5MAXV(20)
      REAL XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,ASM3,VUMT,VDMT,ASMTP,ASMSS,M3Q,MHDSQ,MHUSQ,
     $MHDSMG,MHUSMG,MUMG,BMG,FT2Z1,FB2Z1,FL2Z1,
     $SIGDMX,SIGUMX,C5MAX,C5MAXV
      INTEGER NOGOOD,IAL3UN,ITACHY,MHPNEG,MHLNEG,MHCNEG,
     $MSQNEG,IGUTST
      SAVE /SUGPAS/
C     XNUSUG contains non-universal GUT scale soft terms for SUGRA:
C     XNUSUG(1)=M1 XNUSUG(2)=M2 XNUSUG(3)=M3
C     XNUSUG(4)=A_tau XNUSUG(5)=A_b XNUSUG(6)=A_t
C     XNUSUG(7)=m_Hd XNUSUG(8)=m_Hu XNUSUG(9)=m_eR XNUSUG(10)=m_eL
C     XNUSUG(11)=m_dR XNUSUG(12)=m_uR XNUSUG(13)=m_uL XNUSUG(14)=m_lR
C     XNUSUG(15)=m_lL XNUSUG(16)=m_bR XNUSUG(17)=m_tR XNUSUG(18)=m_tL
C     XNUSUG(19)=mu(Q) XNUSUG(20)=mA(Q)
      COMMON /SUGNU/ XNUSUG(20),INUHM
      REAL XNUSUG
      INTEGER INUHM
      SAVE /SUGNU/
      COMMON /BSG/ GISA(31),MSQISA(3),MSLISA(3),MSUISA(3),MSDISA(3),
     &            MSEISA(3),MRNISA(3),YNFRZ(3,3),MNFRZ(3,3),
     &            TNFRZ(3,3),RTISA,RBISA,RLISA
c     MSxDEC(i) - decoupling scale of i-th generation of type x sfermion
c     MRNDEC(i) - decoupling scale of i-th RH neutrino
      REAL*8 GISA,MSQISA,MSLISA,MSUISA,MSDISA,MSEISA,MRNISA,
     &       YNFRZ,MNFRZ,TNFRZ
      REAL RTISA,RBISA,RLISA
      SAVE /BSG/
      REAL*8 GY(9),W1(27),G(31),W2(93),DT,T,DPI,DY
      REAL*8 BTHAT,BBHAT,BLHAT
      REAL G0(31)
      COMPLEX*16 SSB0,SSB1
      DOUBLE PRECISION DDILOG,XLM
      INTEGER IG(31)
      EXTERNAL SURG06,SURG26
      REAL M0,MHF,A0,TANB,SGNMU,MT,XLAMGM,XMESGM,XN5GM
      INTEGER NSTEP
      REAL M2,SUALFE,SUALFS,Q,A1I,AGUT,A3I,A2I,MTMT,ASMT,
     $TGUT,TZ,GGUT,SIG2,SIG1,SIG3,MH1S,MH2S,AGUTI,
     $MUS,MBMZ,MB,MTAU,MZ,MW,SR2,PI,ALEM,MTAMZ,TANBQ,SIN2BQ,
     $MTAMB,MTAMTA,MBMB,ASMB,BETA,COTB,SINB,COS2B,COSB,XC,
     $L1,L2,L3
      REAL XTGSS,ATGSS,MUSGSS,MUGSS,MGLGSS,MBMZC,MZQ,SIGA
      INTEGER II,I,J,IMODEL
      REAL G0SAVE(31),DELG0,DEL,DELLIM(31),THRF,THRG,QNEW
      REAL CF,CA,ZETA2,ZETA3,ST2LP
      INTEGER MXITER,NSTEP0,IG0LIM
      LOGICAL BADMU
C
      DATA MZ/91.187/,MTAU/1.777/,MB/4.9/,ALEM/.0078186/
      DATA ZETA2/1.644934/,ZETA3/1.202057/
C          This choice is a compromise between precision and speed:
      DATA MXITER/25/,NSTEP0/1000/
C          The Higgs-related soft couplings converge much more slowly
C          than the others, so we use different error cuts:
      DATA DELLIM/
     $0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003,
     $0.003, 0.003, 0.003, 0.003, 0.030, 0.030, 0.003, 0.003,
     $0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003,
     $0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050/
C
C          Define REAL(COMPLEX*16) for g77. This might need to be
C          changed for 64-bit machines?
C
C          Save input parameters
C
      XSUGIN(1)=M0
      XSUGIN(2)=MHF
      XSUGIN(3)=A0
      XSUGIN(4)=TANB
      XSUGIN(5)=SGNMU
      XSUGIN(6)=MT
      XLAMGM=M0
      XMESGM=MHF
      XN5GM=A0
      XGMIN(1)=XLAMGM
      XGMIN(2)=XMESGM
      XGMIN(3)=XN5GM
      XGMIN(4)=TANB
      XGMIN(5)=SGNMU
      XGMIN(6)=MT
      IF (INUHM.EQ.1) THEN
        MU=XNUSUG(19)
        AMHA=XNUSUG(20)
      END IF
      IF (XGMIN(12).EQ.0.) XGMIN(12)=XN5GM
      IF (XGMIN(13).EQ.0.) XGMIN(13)=XN5GM
      IF (XGMIN(14).EQ.0.) XGMIN(14)=XN5GM
      GORGE=.TRUE.
C
C          Compute gauge mediated threshold functions
C
      IF (IMODEL.EQ.2) THEN
        XLM=XLAMGM/XMESGM
        THRF=((1.D0+XLM)*(LOG(1.D0+XLM)-2*DDILOG(XLM/(1.D0+XLM))+
     ,        .5*DDILOG(2*XLM/(1.D0+XLM)))+
     ,       (1.D0-XLM)*(LOG(1.D0-XLM)-2*DDILOG(-XLM/(1.D0-XLM))+
     ,        .5*DDILOG(-2*XLM/(1.D0-XLM))))/XLM**2
        THRG=((1.D0+XLM)*LOG(1.D0+XLM)+(1.D0-XLM)*LOG(1.D0-XLM))
     , /XLM**2
      END IF
C
C          Initialize standard model parameters in /SSSM/:
C
      AMUP=0.0056
      AMDN=0.0099
      AMST=0.199
      AMCH=1.35
      AMBT=4.9
      AMTP=MT
      AMT=MT
      AME=0.511E-3
      AMMU=0.105
      AMTAU=1.777
      AMZ=91.17
      GAMW=2.12
      GAMZ=2.487
      ALFAEM=1./128.
      SN2THW=0.231
      ALFA2=ALFAEM/SN2THW
      ALQCD4=0.177
      ALFA3=0.118
C
      NOGOOD=0
      ITACHY=0
      IGUTST=0
      PI=4.*ATAN(1.)
      DPI=4.D0*DATAN(1.D0)
      CF=4./3.
      CA=3.
      SR2=SQRT(2.)
C      XW=.2324-1.03E-7*(MT**2-138.**2)
      XW=.231
      MW=MZ*SQRT(1.-XW)
      AMW=MW
      A1MZ=5*ALEM/3./(1.-XW)
      A2MZ=ALEM/XW
      G2=SQRT(4*PI*A2MZ)
      GP=SQRT(3./5.*A1MZ*4.*PI)
      XTANB=TANB
      COTB=1./TANB
      BETA=ATAN(TANB)
      SINB=SIN(BETA)
      COSB=COS(BETA)
      SIN2B=SIN(2*BETA)
      COS2B=COS(2*BETA)
      IF (IMODEL.EQ.1) THEN
        MSUSY=SQRT(M0**2+4*MHF**2)
      ELSE IF (IMODEL.EQ.2) THEN
        MSUSY=XLAMGM/100.
      ELSE IF (IMODEL.EQ.7.OR.IMODEL.EQ.9.OR.IMODEL.EQ.10) THEN
        MSUSY=SQRT(M0**2+(.01*MHF)**2)
      END IF
C     USE PIERCE PRESCRIPTION FOR MAGNITUDE OF VEV
      VEV=(248.6+0.9*LOG(MSUSY/AMZ))/SR2
      V=SQRT(VEV**2/(1.+COTB**2))
C     PREVIOUS PRESCRIPTION
C      V=SQRT(2*MW**2/G2**2/(1.+COTB**2))
      VP=V/TANB
C      VEV=SQRT(V**2+VP**2)
C
C          Compute m(tau), m(b) at z scale using qcd, qed
C          Update to DRbar masses used by Pierce et al.
C
C      MTAMTA=MTAU*(1.-SUALFE(MTAU**2)/PI)
C      MTAMB=MTAMTA*(SUALFE(MB**2)/SUALFE(MTAU**2))**(-27./76.)
C      MTAMZ=MTAMB*(SUALFE(MZ**2)/SUALFE(MB**2))**(-27./80.)
      MTAMZ=1.7463
      FTAMZ=MTAMZ/COSB/VEV
C      ASMB=SUALFS(MB**2,.36,MT,3)
C      MBMB=MB*(1.-4*ASMB/3./PI)
C      ASMZ=SUALFS(MZ**2,.36,MT,3)
C      MBMZ=MBMB*(ASMZ/ASMB)**(12./23.)*
C     $      (SUALFE(MZ**2)/SUALFE(MB**2))**(-3./80.)
      MBMZ=2.83
      ASMT=SUALFS(MT**2,.36,MT,3)
      ST2LP=CF*(ASMT/4./PI)**2*(-43.-12*ZETA2+CF*(-59./8.+30*ZETA2-
     ,48*LOG(2.)*ZETA2+12*ZETA3)+
     ,CA*(1093./24.-8*ZETA2+24*LOG(2.)*ZETA2-6*ZETA3))
      MTMT=MT/(1.+5*ASMT/3./PI+ST2LP)
      FTMT=MTMT/SINB/VEV
C     Here we input Drees' guesses for Z-scale soft terms
C     so we have a good initial guess for mb(mz)
C     Guesses come from hep-ph/9504324
C      IF (IMODEL.EQ.1) THEN
C      XTGSS=(MTMT/150./SINB)**2*(.9*M0**2+2.1*MHF**2+
C     $(1.-(MTMT/190./SINB)**3)*(.24*A0**2+A0*MHF))
C      ATGSS=A0*(1.-(MTMT/190./SINB)**2)-
C     $MHF*(3.47-1.9*(MTMT/190./SINB)**2)
C      MUSGSS=-M0**2-.52*MHF**2-.5*MZ**2+XTGSS/(1.-COTB**2)
C      MUGSS=SQRT(MAX(0.,MUSGSS))*SIGN(1.,SGNMU)
C      MGLGSS=MHF*ASMT/0.04
C      MBMZC=MBMZ*(1.+2*ASMZ/3./PI*MUGSS*MGLGSS/MSUSY**2*TANB+
C     $FTMT**2/16./PI**2*MUGSS*ATGSS/MSUSY**2)
C      ELSE
C      MBMZC=MBMZ
C      END IF
C      FBMZ=MBMZC/COSB/VEV
      FBMZ=MBMZ/COSB/VEV
      FNMZ=SQRT(XNRIN(2)*XNRIN(1)/(SINB*VEV)**2)
      AMNRMJ=XNRIN(2)
C     Initialize some parameters for SUGFRZ
      IF (INUHM.NE.1) AMHA=AMZ
      ASM3=ALFA3
C     Set GSS values to initial guess
      DO I=7,12
        GSS(I)=MSUSY
      END DO
      DO I=13,24
        GSS(I)=MSUSY**2
      END DO
C
C          Run the 3 gauge and 3 Yukawa's up to find M_GUT ,A_GUT and
C          Yukawa_GUT
C
      NSTEP=NSTEP0
      GY(1)=DBLE(SQRT(4*PI*A1MZ))
      GY(2)=DBLE(SQRT(4*PI*A2MZ))
      GY(3)=DBLE(SQRT(4*PI*ALFA3))
      GY(4)=DBLE(FTAMZ)
      GY(5)=DBLE(FBMZ)
      GY(6)=0.D0
      GY(7)=0.D0
      GY(8)=DBLE(VP)
      GY(9)=DBLE(V)
      IF (IMODEL.EQ.1.OR.IMODEL.EQ.7.OR.IMODEL.EQ.9.OR.
     ,IMODEL.EQ.10) THEN
        IF (XSUGIN(7).EQ.0.) THEN
          MGUT=1.E19
        ELSE
          MGUT=XSUGIN(7)
        END IF
      ELSE IF (IMODEL.EQ.2) THEN
        MGUT=XMESGM
      END IF
      TZ=DLOG(DBLE(MZ)/DBLE(MGUT))
      TGUT=0.D0
      DT=(TGUT-TZ)/DBLE(FLOAT(NSTEP))
      DO 200 II=1,NSTEP
        T=TZ+(TGUT-TZ)*FLOAT(II-1)/FLOAT(NSTEP)
        Q=MGUT*EXP(SNGL(T))
        IF (Q.GT.MT.AND.GY(6).EQ.0.D0) GY(6)=DBLE(FTMT)
        IF (Q.GT.XNRIN(2).AND.GY(7).EQ.0.D0) GY(7)=DBLE(FNMZ)
        CALL DRKSTP(9,DT,T,GY,SURG06,W1)
        A1I=4*PI/SNGL(GY(1)**2)
        A2I=4*PI/SNGL(GY(2)**2)
        A3I=4*PI/SNGL(GY(3)**2)
        IF (GY(4).GT.5.D0.OR.GY(5).GT.5.D0.OR.
     $  GY(6).GT.5.D0.OR.GY(7).GT.5.D0) THEN
          NOGOOD=4
          GO TO 100
        END IF
        IF (A1I.LT.A2I.AND.XSUGIN(7).EQ.0.) GO TO 10
200   CONTINUE
      IF (MGUT.EQ.1.E19) THEN
        WRITE(LOUT,*) 'SUGRA: NO UNIFICATION FOUND'
        GO TO 100
      END IF
10    IF (XSUGIN(7).EQ.0.) THEN
        MGUT=Q
      ELSE
        MGUT=XSUGIN(7)
      END IF
      AGUT=SNGL((GY(1)**2/4.D0/DPI+GY(2)**2/4.D0/DPI)/2.D0)
      GGUT=SQRT(4*PI*AGUT)
      AGUTI=1./AGUT
      FTAGUT=SNGL(GY(4))
      FBGUT=SNGL(GY(5))
      FTGUT=SNGL(GY(6))
      IF (XNRIN(1).EQ.0..AND.XNRIN(2).LT.1.E19) THEN
C       UNIFY FN-FT
        FNGUT=SNGL(GY(6))
      ELSE
        FNGUT=SNGL(GY(7))
      END IF
C
C          Define parameters at GUT scale
C
      DO 210 J=1,3
        IF (IMODEL.EQ.1) THEN
          G(J)=GY(J)
          G(J+6)=DBLE(MHF)
          G(J+9)=DBLE(A0)
        ELSE IF (IMODEL.EQ.2) THEN
          G(J)=GY(J)
          G(J+6)=DBLE(XGMIN(11+J)*XGMIN(8)*THRG*(GY(J)/4./PI)**2
     &*XLAMGM)
          G(J+9)=0.D0
        END IF
210   CONTINUE
      G(30)=GY(8)
      G(31)=GY(9)
C          Overwrite alfa_3 unification to get alfa_3(mz) right
      IF (IMODEL.EQ.1.AND.IAL3UN.NE.0) G(3)=DBLE(GGUT)
      G(4)=DBLE(FTAGUT)
      G(5)=DBLE(FBGUT)
      G(6)=DBLE(FTGUT)
C          If nr Majorana mass exists, set extra nr rge parameters
      IF (XNRIN(2).LT.1.E19) THEN
        G(27)=DBLE(FNGUT)
        G(28)=DBLE(XNRIN(4))**2
        G(29)=DBLE(XNRIN(3))
      ELSE
        G(27)=0.D0
        G(28)=0.D0
        G(29)=0.D0
      END IF
      IF (IMODEL.EQ.1) THEN
        DO 220 J=13,24
          G(J)=DBLE(M0)**2
220     CONTINUE
C     Azar suggested NUHM2 boundary condition guess
        IF (INUHM.EQ.1) THEN
c...compute approximate MSUSY-scale higgs mass parameters
          MHDSQ=(-MZ**2-2*MU**2+2*TANB**2*(AMHA**2+MZ**2)
     &           +TANB**4*(2*MU**2-MZ**2-2*AMHA**2))/2./(1.-TANB**4)
          MHUSQ=(2*AMHA**2+MZ**2-2*MU**2+2*TANB**2*(-AMHA**2-MZ**2)
     &           +TANB**4*(2*MU**2+MZ**2))/2./(1.-TANB**4)
c...use analytical 1lp solutions from hep-ph/0103324 for tanb=10
c   to find approximate GUT-scale values
        MHUSMG=(MHUSQ+0.75*M0**2+2.71*MHF**2-0.39*A0*MHF)/0.63
        MHDSMG=(MHDSQ-0.21*MHF**2+0.02*M0**2)/0.99
c...add dominant 2-loop piece to up Higgs
        if(XNUSUG(13).LT.1.E19) THEN
         MHUSMG=MHUSMG
     &       +2*12./5./(16.*PI**2)**2*(G(1)**4+9.*G(2)**2)
     &*XNUSUG(13)**2
     &        *log(XNUSUG(13)/MGUT)
          else
         MHUSMG=MHUSMG
     &          +3*12./5./(16.*PI**2)**2*(G(1)**4+9.*G(2)**2)*m0**2
     &           *log(m0/MGUT)
        endif
c        print*,'MHUSMG,MHDSMG=',MHUSMG,MHDSMG
c          G(13)=DBLE(MHDSMG)
c          G(14)=DBLE(MHUSMG)
         ENDIF
C       End Azar guess
C
C       Set possible non-universal boundary conditions
        DO 230 J=1,6
          IF (XNUSUG(J).LT.1.E19) THEN
            G(J+6)=DBLE(XNUSUG(J))
          END IF
230     CONTINUE
        DO 231 J=7,18
          IF (XNUSUG(J).LT.1.E19) THEN
            G(J+6)=SIGN(1.,XNUSUG(J))*DBLE(XNUSUG(J))**2
          END IF
231     CONTINUE
      ELSE IF (IMODEL.EQ.2) THEN
        XC=2*THRF*XLAMGM**2
        DY=DSQRT(3.D0/5.D0)*GY(1)*XGMIN(11)
        G(13)=XC*(.75*XGMIN(13)*(GY(2)/4.D0/DPI)**4+.6D0*.25*
     ,  XGMIN(12)*(GY(1)/4.D0/DPI)**4)+DBLE(XGMIN(9))-DY
        G(14)=XC*(.75*XGMIN(13)*(GY(2)/4.D0/DPI)**4+.6D0*.25*
     ,  XGMIN(12)*(GY(1)/4.D0/DPI)**4)+DBLE(XGMIN(10))+DY
        G(15)=XC*(.6*XGMIN(12)*(GY(1)/4.D0/DPI)**4)+2*DY
        G(16)=XC*(.75*XGMIN(13)*(GY(2)/4.D0/DPI)**4+.6D0*.25*
     ,  XGMIN(12)*(GY(1)/4.D0/DPI)**4)-DY
        G(17)=XC*(4*XGMIN(14)*(GY(3)/4.D0/DPI)**4/3.D0+.6D0*
     ,  XGMIN(12)*(GY(1)/4.D0/DPI)**4/9.D0)+2*DY/3.D0
        G(18)=XC*(4*XGMIN(14)*(GY(3)/4.D0/DPI)**4/3.D0+.6D0*
     ,  4*XGMIN(12)*(GY(1)/4.D0/DPI)**4/9.D0)-4*DY/3.D0
        G(19)=XC*(4*XGMIN(14)*(GY(3)/4.D0/DPI)**4/3.D0+.75D0*
     ,  XGMIN(13)*(GY(2)/4.D0/DPI)**4+.6D0*XGMIN(12)*(GY(1)/
     ,  4.D0/DPI)**4/36.D0)+DY/3.D0
        G(20)=G(15)
        G(21)=G(16)
        G(22)=G(17)
        G(23)=G(18)
        G(24)=G(19)
      ELSE IF (IMODEL.EQ.7.OR.IMODEL.EQ.9.OR.IMODEL.EQ.10) THEN
        G(1)=GY(1)
        G(2)=GY(2)
        G(3)=GY(3)
        BLHAT=G(4)*(-9*G(1)**2/5.D0-3*G(2)**2+3*G(5)**2+4*G(4)**2)
        BBHAT=G(5)*(-7*G(1)**2/15.D0-3*G(2)**2-16*G(3)**2/3.D0+
     ,             G(6)**2+6*G(5)**2+G(4)**2)
        BTHAT=G(6)*(-13*G(1)**2/15.D0-3*G(2)**2-16*G(3)**2/3.D0+
     ,             6*G(6)**2+G(5)**2)
        G(7)=33*MHF*G(1)**2/5.D0/16.D0/DPI**2
        IF (IMODEL.EQ.10) THEN
          G(7)=G(7)+XAMIN(11)*MHF
        END IF
        G(8)=MHF*G(2)**2/16.D0/DPI**2
        G(9)=-3*MHF*G(3)**2/16.D0/DPI**2
        G(10)=-BLHAT*MHF/G(4)/16.D0/DPI**2
        G(11)=-BBHAT*MHF/G(5)/16.D0/DPI**2
        G(12)=-BTHAT*MHF/G(6)/16.D0/DPI**2
        G(13)=(-99*G(1)**4/50.D0-3*G(2)**4/2.D0+3*G(5)*BBHAT+
     ,G(4)*BLHAT)*MHF**2/(16*DPI**2)**2+XAMIN(6)*DBLE(M0)**2
        G(14)=(-99*G(1)**4/50.D0-3*G(2)**4/2.D0+3*G(6)*BTHAT)*
     ,        MHF**2/(16*DPI**2)**2+XAMIN(7)*DBLE(M0)**2
        G(15)=(-198*G(1)**4/25.D0)*MHF**2/(16*DPI**2)**2+
     ,XAMIN(5)*DBLE(M0)**2
        G(16)=(-99*G(1)**4/50.D0-3*G(2)**4/2.D0)*MHF**2/(16*DPI**2)**2
     ,+XAMIN(4)*DBLE(M0)**2
        G(17)=(-22*G(1)**4/25.D0+8*G(3)**4)*MHF**2/(16*DPI**2)**2+
     ,XAMIN(2)*DBLE(M0)**2
        G(18)=(-88*G(1)**4/25.D0+8*G(3)**4)*MHF**2/(16*DPI**2)**2+
     ,XAMIN(3)*DBLE(M0)**2
        G(19)=(-11*G(1)**4/50.D0-3*G(2)**4/2.D0+8*G(3)**4)*
     ,        MHF**2/(16*DPI**2)**2+XAMIN(1)*DBLE(M0)**2
        G(20)=(-198*G(1)**4/25.D0+2*G(4)*BLHAT)*MHF**2/(16*DPI**2)**2
     ,+XAMIN(5)*DBLE(M0)**2
        G(21)=(-99*G(1)**4/50.D0-3*G(2)**4/2.D0+G(4)*BLHAT)*
     ,        MHF**2/(16*DPI**2)**2+XAMIN(4)*DBLE(M0)**2
        G(22)=(-22*G(1)**4/25.D0+8*G(3)**4+2*G(5)*BBHAT)*
     ,        MHF**2/(16*DPI**2)**2+XAMIN(2)*DBLE(M0)**2
        G(23)=(-88*G(1)**4/25.D0+8*G(3)**4+2*G(6)*BTHAT)*
     ,        MHF**2/(16*DPI**2)**2+XAMIN(3)*DBLE(M0)**2
        G(24)=(-11*G(1)**4/50.D0-3*G(2)**4/2.D0+8*G(3)**4+G(5)*BBHAT+
     ,        G(6)*BTHAT)*MHF**2/(16*DPI**2)**2+XAMIN(1)*DBLE(M0)**2
      END IF
      IF (IMODEL.EQ.9) THEN
        CALL MMAMSB(M0,MHF,G)
      END IF
      G(25)=0.D0
      G(26)=0.D0
      DO 235 I=1,31
        IG(I)=0
235   CONTINUE
C          Check for tachyonic sleptons at GUT scale
      IF (G(15).LT.0.D0.OR.G(16).LT.0.D0) THEN
        ITACHY=1
      END IF
C
C          Initialize all masses to MSUSY scale
C
      DO 236 I=1,32
        MSS(I)=MSUSY+FLOAT(I)
236   CONTINUE
      IF (INUHM.EQ.1) THEN
        MSS(31)=AMHA
      END IF
C
C          Evolve parameters from mgut to mz
C
      TZ=DLOG(DBLE(MZ)/DBLE(MGUT))
      TGUT=0.D0
      DT=(TZ-TGUT)/DBLE(FLOAT(NSTEP))
C          Freeze Higgs parameters at HIGFRZ = Drees' value
C          AMTLSS, AMTRSS initialized to 0 for later use in HIGFRZ
      IF (IMODEL.EQ.1) THEN
        HIGFRZ=SQRT(M0**2+3*MHF**2)
      ELSE IF (IMODEL.EQ.2) THEN
        HIGFRZ=MSUSY
      ELSE IF (IMODEL.EQ.7.OR.IMODEL.EQ.9.OR.IMODEL.EQ.10) THEN
        HIGFRZ=SQRT(M0**2+(.01*MHF)**2)
      END IF
      AMTLSS=0.
      AMTRSS=0.
      DO 240 II=1,NSTEP+2
        T=TGUT+(TZ-TGUT)*FLOAT(II-1)/DBLE(FLOAT(NSTEP))
        Q=SNGL(MGUT*DEXP(T))
        CALL DRKSTP(31,DT,T,G,SURG26,W2)
        QNEW=SNGL(MGUT*DEXP(T+DT))
C       TEST YUKAWA DIVERGENCE
        IF (G(4).GT.5.D0.OR.G(5).GT.5.D0.OR.
     $  G(6).GT.5.D0.OR.G(27).GT.5.D0) THEN
          NOGOOD=4
          GO TO 100
        END IF
        IF (QNEW.LT.AMNRMJ.AND.Q.GE.AMNRMJ.AND.FNMZ.EQ.0.) THEN
          FNMZ=SNGL(G(27))
        END IF
        IF (QNEW.LT.AMNRMJ) THEN
          G(27)=0.D0
          G(28)=0.D0
          G(29)=0.D0
        END IF
        CALL SUGFRZ(QNEW,G,G0,IG)
        IF (NOGOOD.NE.0) GO TO 100
        IF (QNEW.LT.MZ) GO TO 20
240   CONTINUE
20    CONTINUE
      ASMZ=G0(3)**2/4./PI
      VUQ=G0(31)
      VDQ=G0(30)
      TANBQ=VUQ/VDQ
      SIN2BQ=SIN(2*ATAN(TANBQ))
      MZQ=SQRT((G0(2)**2+.6*G0(1)**2)*(VUQ**2+VDQ**2)/2.)
      IF (INUHM.NE.1) THEN
C          Electroweak breaking constraints
C          Tree level
      MUS=(G0(13)-G0(14)*TANBQ**2)/(TANBQ**2-1.)-MZQ**2/2.
C          Calculate loop corrections using MSS=MSUSY masses set above
      CALL SUGEFF(G0,SIG1,SIG2,SIG3)
      MH1S=G0(13)+SIG1
      MH2S=G0(14)+SIG2
      MUS=(MH1S-MH2S*TANBQ**2)/(TANBQ**2-1.)-MZQ**2/2.
C          If MUS<0, set it to MZ**2 so that spectra and real loop
C          corrections can be calculated
      IF (MUS.LT.0.) THEN
        MUS=MZ**2
      END IF
      MU=SQRT(MUS)*SIGN(1.,SGNMU)
      B=(G0(13)+G0(14)+2*MUS)*SIN2BQ/MU/2.+SIG3/MU
C          Compute tree level masses using first value of MU
      CALL SUGMAS(G0,0,IMODEL,SIGA)
      IF (NOGOOD.NE.0) GO TO 100
C          Compute effective potential corrections with tree masses
      CALL SUGEFF(G0,SIG1,SIG2,SIG3)
      MH1S=G0(13)+SIG1
      MH2S=G0(14)+SIG2
      MUS=(MH1S-MH2S*TANBQ**2)/(TANBQ**2-1.)-MZQ**2/2.
C          MUS might still be negative. If so, set it to MZ**2 and hope
C          for the best....
      IF (MUS.LT.0.) THEN
C       NOGOOD=2
C       GO TO 100
        MUS=MZ**2
      END IF
      MU=SQRT(MUS)*SIGN(1.,SGNMU)
      B=(MH1S+MH2S+2*MUS)*SIN2BQ/MU/2.+SIG3/MU
C          Need loop corrected mass spectra to calculate fermion
C          self energies
      CALL SUGMAS(G0,1,IMODEL,SIGA)
      ELSE
      MUS=MU**2
      B=AMHA**2/MU/(COTB+TANB)+SIG3/MU
      CALL SUGMAS(G0,1,IMODEL,SIGA)
      CALL SUGEFF(G0,SIG1,SIG2,SIG3)
      MHDSQ=(AMHA**2+MZQ**2)*TANB**2/(TANB**2+1.)-SIG1-MUS-MZQ**2/2.
      MHUSQ=(AMHA**2+MZQ**2)/(TANB**2+1.)-MUS-MZQ**2/2.-SIG2
      END IF
C
C          Iterate entire process, increasing NSTEP each time
C          This time, freeze out parameters at sqrt(t_l t_r)
C
      HIGFRZ=(MAX(AMZ**4,G0(23)*G0(24)))**0.25
      MSUSY=HIGFRZ
      DO 300 I=1,MXITER
        DO 310 J=1,31
310     G0SAVE(J)=G0(J)
        NSTEP=1.2*NSTEP
c        CALL SUGRGE(M0,MHF,A0,TANB,SGNMU,MT,G,G0,IG,W2,NSTEP
c     $ ,IMODEL,BADMU)
        CALL MODSUGRGE(M0,MHF,A0,TANB,SGNMU,MT,G,G0,IG,W2,NSTEP
     $ ,IMODEL,BADMU)
        HIGFRZ=(MAX(AMZ**4,G0(23)*G0(24)))**0.25
        MSUSY=HIGFRZ
        IF(NOGOOD.NE.0) GO TO 100
C            Check convergence relative to DELLIM
        DELG0=0.
        IG0LIM=0
        DO 320 J=1,31
          IF(G0(J).NE.0) THEN
            DEL=ABS((G0(J)-G0SAVE(J))/G0(J))
          ELSE
            DEL=0
          ENDIF
          IF(DEL-DELLIM(J).GT.DELG0) THEN
            DELG0=DEL
            IG0LIM=J
          ENDIF
320     CONTINUE
C       Azar's GSS fix
        DO J=1,31
          GSS(J)=G0(J)
        ENDDO
        IF(IG0LIM.EQ.0) GO TO 400
300   CONTINUE
C
C          No solution found in MXITER iterations
C
      WRITE(LOUT,1000) MXITER,DELG0,IG0LIM
1000  FORMAT(/' SUGRA: NO RGE CONVERGENCE IN',I4,' ITERATIONS'/
     $' WORST ERROR = ',E12.4,' FOR G0(',I2,')')
      NOGOOD=-1
      GO TO 100
C
C          Save results
C
400   DO 410 I=1,31
        GSS(I)=G0(I)
410   CONTINUE
C          Set flag for NOGOOD radiative EWSB. We allow MUS to be
C          negative at intermediate steps but check here that it is
C          positive after all iterations.
      IF (BADMU) THEN
        NOGOOD=2
      END IF
C          Set flag for NOGOOD GUT stability in NUHM model
      IF (INUHM.EQ.1) THEN
      IF (MHUSMG+MUMG**2.LT.0..OR.MHDSMG+MUMG**2.LT.0.) THEN
        IGUTST=1
      END IF
      END IF
      MGUTSS=MGUT
      AGUTSS=AGUT
      GGUTSS=GGUT
C     Compute nu_3 mass
      IF (XNRIN(1).EQ.0..AND.XNRIN(2).LT.1.E19) THEN
         XNRIN(1)=(GSS(27)*SINB*VEV)**2/GSS(28)
      ENDIF
C
C          Fill XISAIN common block
C
      XISAIN(1)=MSS(1)
      XISAIN(2)=MU
      XISAIN(3)=MSS(31)
      XISAIN(4)=TANB
      XISAIN(5)=SQRT(G0(19))
      XISAIN(6)=SQRT(G0(17))
      XISAIN(7)=SQRT(G0(18))
      XISAIN(8)=SQRT(G0(16))
      XISAIN(9)=SQRT(G0(15))
      XISAIN(10)=XISAIN(5)
      XISAIN(11)=XISAIN(6)
      XISAIN(12)=XISAIN(7)
      XISAIN(13)=XISAIN(8)
      XISAIN(14)=XISAIN(9)
C     KEEP TRACK OF SIGN OF SQUARED SSB TERMS MTL**2 AND MTR**2
      XISAIN(15)=SIGN(1.,G0(24))*SQRT(ABS(G0(24)))
      XISAIN(16)=SQRT(G0(22))
      XISAIN(17)=SIGN(1.,G0(23))*SQRT(ABS(G0(23)))
      XISAIN(18)=SQRT(G0(21))
      XISAIN(19)=SQRT(G0(20))
      XISAIN(20)=G0(12)
      XISAIN(21)=G0(11)
      XISAIN(22)=G0(10)
      XISAIN(23)=G0(7)
      XISAIN(24)=G0(8)
      M2=G0(8)
c
c   Save values of RGE parameters at MZ for b->s\gamma computation
c
      DO 420 I=1,31
        GISA(I)=G(I)
420   CONTINUE
      MSQISA(1)=DBLE(MSS(2))
      MSQISA(2)=DBLE(MSS(2))
      MSQISA(3)=DBLE(MSS(2))
      MSLISA(1)=DBLE(MSS(17))
      MSLISA(2)=DBLE(MSS(17))
      MSLISA(3)=DBLE(MSS(17))
      MSUISA(1)=DBLE(MSS(2))
      MSUISA(2)=DBLE(MSS(2))
      MSUISA(3)=DBLE(MSS(2))
      MSDISA(1)=DBLE(MSS(2))
      MSDISA(2)=DBLE(MSS(2))
      MSDISA(3)=DBLE(MSS(2))
      MSEISA(1)=DBLE(MSS(17))
      MSEISA(2)=DBLE(MSS(17))
      MSEISA(3)=DBLE(MSS(17))
      MRNISA(3)=AMNRMJ
      MRNISA(1)=MRNISA(3)
      MRNISA(2)=MRNISA(3)
C
100   RETURN
      END



C....Modified version of subroutine SUGRGE. Save running parameters in block TRHPARS
      SUBROUTINE MODSUGRGE(M0,MHF,A0,TANB,SGNMU,MT,G,G0,IG,W2
     $,NSTEP,IMODEL,BADMU)
C
C          Make one complete iteration of the renormalization group
C          equations from MZ to MGUT and back, setting the boundary
C          conditions on each end.
C
      IMPLICIT NONE


      DOUBLE PRECISION TRHM1V,TRHM2V,TRHM3V,TRHG1V,TRHG2V,TRHG3V,TQV
      COMMON/TRHPARS/ TRHM1V(100000),TRHM2V(100000),TRHM3V(100000),
     &TRHG1V(100000),TRHG2V(100000),TRHG3V(100000),TQV(100000)
      SAVE/TRHPARS/



      COMMON/SSLUN/LOUT,LHEOUT
      INTEGER LOUT,LHEOUT
      SAVE /SSLUN/
C          Standard model parameters
C          AMUP,...,AMTP        = quark masses
C          AME,AMMU,AMTAU       = lepton masses
C          AMW,AMZ              = W,Z masses
C          GAMW,GAMZ            = W,Z widths
C          ALFAEM,SN2THW,ALFA3  = SM couplings
C          ALQCD4               = 4 flavor lambda
      COMMON/SSSM/AMUP,AMDN,AMST,AMCH,AMBT,AMTP,AME,AMMU,AMTAU
     $,AMW,AMZ,GAMW,GAMZ,ALFAEM,SN2THW,ALFA2,ALFA3,ALQCD4
      REAL AMUP,AMDN,AMST,AMCH,AMBT,AMTP,AME,AMMU,AMTAU
     $,AMW,AMZ,GAMW,GAMZ,ALFAEM,SN2THW,ALFA2,ALFA3,ALQCD4
      SAVE /SSSM/
C          SUSY parameters
C          AMGLSS               = gluino mass
C          AMULSS               = up-left squark mass
C          AMELSS               = left-selectron mass
C          AMERSS               = right-slepton mass
C          AMNiSS               = sneutrino mass for generation i
C          TWOM1                = Higgsino mass = - mu
C          RV2V1                = ratio v2/v1 of vev's
C          AMTLSS,AMTRSS        = left,right stop masses
C          AMT1SS,AMT2SS        = light,heavy stop masses
C          AMBLSS,AMBRSS        = left,right sbottom masses
C          AMB1SS,AMB2SS        = light,heavy sbottom masses
C          AMLLSS,AMLRSS        = left,right stau masses
C          AML1SS,AML2SS        = light,heavy stau masses
C          AMZiSS               = signed mass of Zi
C          ZMIXSS               = Zi mixing matrix
C          AMWiSS               = signed Wi mass
C          GAMMAL,GAMMAR        = Wi left, right mixing angles
C          AMHL,AMHH,AMHA       = neutral Higgs h0, H0, A0 masses
C          AMHC                 = charged Higgs H+ mass
C          ALFAH                = Higgs mixing angle
C          AAT                  = stop trilinear term
C          THETAT               = stop mixing angle
C          AAB                  = sbottom trilinear term
C          THETAB               = sbottom mixing angle
C          AAL                  = stau trilinear term
C          THETAL               = stau mixing angle
C          AMGVSS               = gravitino mass
C          MTQ                  = top mass at MSUSY
C          MBQ                  = bottom mass at MSUSY
C          MLQ                  = tau mass at MSUSY
C          FBMA                 = b-Yukawa at mA scale
C          VUQ                  = Hu vev at MSUSY
C          VDQ                  = Hd vev at MSUSY
C          SGNM3                = sign of gaugino mass M3
      COMMON/SSPAR/GORGE,AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS(4,4)
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,
     $VUQ,VDQ,SGNM3
      REAL AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,VUQ,VDQ,SGNM3
      LOGICAL GORGE
      REAL AMZISS(4)
      EQUIVALENCE (AMZISS(1),AMZ1SS)
      SAVE /SSPAR/
      COMMON /SUGPAS/ XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,NOGOOD,IAL3UN,ITACHY,MHPNEG,MHLNEG,MHCNEG,
     $MSQNEG,IGUTST,ASM3,
     $VUMT,VDMT,ASMTP,ASMSS,M3Q,MHDSQ,MHUSQ,MHDSMG,MHUSMG,MUMG,BMG,
     $FT2Z1,FB2Z1,FL2Z1,SIGDMX,SIGUMX,C5MAX,C5MAXV(20)
      REAL XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,ASM3,VUMT,VDMT,ASMTP,ASMSS,M3Q,MHDSQ,MHUSQ,
     $MHDSMG,MHUSMG,MUMG,BMG,FT2Z1,FB2Z1,FL2Z1,
     $SIGDMX,SIGUMX,C5MAX,C5MAXV
      INTEGER NOGOOD,IAL3UN,ITACHY,MHPNEG,MHLNEG,MHCNEG,
     $MSQNEG,IGUTST
      SAVE /SUGPAS/
C     XNUSUG contains non-universal GUT scale soft terms for SUGRA:
C     XNUSUG(1)=M1 XNUSUG(2)=M2 XNUSUG(3)=M3
C     XNUSUG(4)=A_tau XNUSUG(5)=A_b XNUSUG(6)=A_t
C     XNUSUG(7)=m_Hd XNUSUG(8)=m_Hu XNUSUG(9)=m_eR XNUSUG(10)=m_eL
C     XNUSUG(11)=m_dR XNUSUG(12)=m_uR XNUSUG(13)=m_uL XNUSUG(14)=m_lR
C     XNUSUG(15)=m_lL XNUSUG(16)=m_bR XNUSUG(17)=m_tR XNUSUG(18)=m_tL
C     XNUSUG(19)=mu(Q) XNUSUG(20)=mA(Q)
      COMMON /SUGNU/ XNUSUG(20),INUHM
      REAL XNUSUG
      INTEGER INUHM
      SAVE /SUGNU/
C     XSUGIN contains the inputs to SUGRA:
C     XSUGIN(1) = M_0        XSUGIN(2) = M_(1/2)  XSUGIN(3) = A_0
C     XSUGIN(4) = tan(beta)  XSUGIN(5) = sgn(mu)  XSUGIN(6) = M_t
C     XSUGIN(7) = SUG BC scale
C     XGMIN(1) = LAM         XGMIN(2)  = M_MES    XGMIN(3)  = XN5
C     XGMIN(4) = tan(beta)   XGMIN(5)  = sgn(mu)  XGMIN(6) = M_t
C     XGMIN(7) = CGRAV       XGMIN(8)  =RSL       XGMIN(9)  = DEL_HD
C     XGMIN(10)  = DEL_HU    XGMIN(11) = DY       XGMIN(12) = N5_1
C     XGMIN(13)  = N5_2      XGMIN(14) = N5_3
C     XNRIN(1) = M_N3        XNRIN(2) = M_MAJ     XNRIN(3) = ANSS
C     XNRIN(4) = M_N3SS
C     XISAIN contains the MSSMi inputs in natural order.
      COMMON /SUGXIN/ XISAIN(24),XSUGIN(7),XGMIN(14),XNRIN(4),
     $XAMIN(11)
      REAL XISAIN,XSUGIN,XGMIN,XNRIN,XAMIN
      SAVE /SUGXIN/
C          Frozen couplings from RG equations:
C     GSS( 1) = g_1        GSS( 2) = g_2        GSS( 3) = g_3
C     GSS( 4) = y_tau      GSS( 5) = y_b        GSS( 6) = y_t
C     GSS( 7) = M_1        GSS( 8) = M_2        GSS( 9) = M_3
C     GSS(10) = A_tau      GSS(11) = A_b        GSS(12) = A_t
C     GSS(13) = M_hd^2     GSS(14) = M_hu^2     GSS(15) = M_er^2
C     GSS(16) = M_el^2     GSS(17) = M_dnr^2    GSS(18) = M_upr^2
C     GSS(19) = M_upl^2    GSS(20) = M_taur^2   GSS(21) = M_taul^2
C     GSS(22) = M_btr^2    GSS(23) = M_tpr^2    GSS(24) = M_tpl^2
C     GSS(25) = mu         GSS(26) = B          GSS(27) = Y_N
C     GSS(28) = M_nr       GSS(29) = A_n        GSS(30) = vdq
C     GSS(31) = vuq
C          Masses:
C     MSS( 1) = glss     MSS( 2) = upl      MSS( 3) = upr
C     MSS( 4) = dnl      MSS( 5) = dnr      MSS( 6) = stl
C     MSS( 7) = str      MSS( 8) = chl      MSS( 9) = chr
C     MSS(10) = b1       MSS(11) = b2       MSS(12) = t1
C     MSS(13) = t2       MSS(14) = nuel     MSS(15) = numl
C     MSS(16) = nutl     MSS(17) = el-      MSS(18) = er-
C     MSS(19) = mul-     MSS(20) = mur-     MSS(21) = tau1
C     MSS(22) = tau2     MSS(23) = z1ss     MSS(24) = z2ss
C     MSS(25) = z3ss     MSS(26) = z4ss     MSS(27) = w1ss
C     MSS(28) = w2ss     MSS(29) = hl0      MSS(30) = hh0
C     MSS(31) = ha0      MSS(32) = h+
C          Unification:
C     MGUTSS  = M_GUT    GGUTSS  = g_GUT    AGUTSS  = alpha_GUT
      COMMON /SUGMG/ MSS(32),GSS(31),MGUTSS,GGUTSS,AGUTSS,FTGUT,
     $FBGUT,FTAGUT,FNGUT
      REAL MSS,GSS,MGUTSS,GGUTSS,AGUTSS,FTGUT,FBGUT,FTAGUT,FNGUT
      SAVE /SUGMG/
      COMMON/SSINF/XLAM
      DOUBLE PRECISION XLAM
C
      COMMON /BSG/GISA(31),MSQISA(3),MSLISA(3),MSUISA(3),MSDISA(3),
     &            MSEISA(3),MRNISA(3),YNFRZ(3,3),MNFRZ(3,3),
     &            TNFRZ(3,3),RTISA,RBISA,RLISA
c     MSxDEC(i) - decoupling scale of i-th generation of type x sfermion
c     MRNDEC(i) - decoupling scale of i-th RH neutrino
      REAL*8 GISA,MSQISA,MSLISA,MSUISA,MSDISA,MSEISA,MRNISA,
     &       YNFRZ,MNFRZ,TNFRZ
      REAL RTISA,RBISA,RLISA
      SAVE /BSG/
C     Common blocks for A. Box RGE code
      COMMON/RGEMS/VEVMH,RGEMS,RGEMU
      DOUBLE COMPLEX VEVMH
      DOUBLE PRECISION RGEMS,RGEMU
      SAVE/RGEMS/
C
      COMMON/WKYUK/LAMTMT,LAMBMZ,LAMTAMZ
      DOUBLE PRECISION LAMTMT,LAMBMZ,LAMTAMZ
      SAVE/WKYUK/
C
C
      EXTERNAL SURG26
      DOUBLE PRECISION DDILOG,XLM
      COMPLEX*16 SSB0,SSB1
      REAL*8 G(31),W2(93),T,DT,DY,DPI,DAS
      REAL*8 BTHAT,BBHAT,BLHAT
      REAL M0,MHF,A0,TANB,SGNMU,MT,G0(31)
      INTEGER IG(31),NSTEP,IMODEL
      REAL PI,TZ,A1I,A2I,A3I,GGUT,AGUTI,SIG1,SIG2,SIG3,
     $MH1S,MH2S,MUS,MZ,TGUT,AGUT,Q,ASMT,MTMT,
     $QNEW,XLAMGM,XMESGM,XN5GM,XC,G3GUT,THRF,THRG,MBMZ,
     $M2,AM2,MSN,MG,MT1,MT2,MB1,MB2,MW1,MW2,AMU,
     $RSIGT,RSIGL,RSIGB,DAEM,ALEMDR,ALEM,MTAMZ,
     $TANBQ,SIN2BQ,SINB,COSB,XWMSB
      REAL SSRSGT,SSRSGB,SSRSGL,SUALFS,MZQ,SIGA
      REAL CF,CA,ZETA2,ZETA3,ST2LP,COTB
      INTEGER I,II
      LOGICAL BADMU
C
      DATA ZETA2/1.644934/,ZETA3/1.202057/
      DATA MZ/91.187/
C
C          Re-initialize weak scale parameters
C
      XLAMGM=M0
      XMESGM=MHF
      XN5GM=A0
      PI=4.*ATAN(1.)
      DPI=4.D0*DATAN(1.D0)
      CF=4./3.
      CA=3.
C     Here we input alpha_s^MSbar(MZ)
      ASMZ=0.1172
      MTAMZ=1.7463
C     Value of mb(MZ)^DRbar taken from PRD66, 074007 (2002).
      MBMZ=2.83
      SINB=SIN(ATAN(TANB))
      COSB=COS(ATAN(TANB))
      COTB=1./TANB
C
C     Calculate fermion masses including loop corrections
C
      M2=G0(8)
      AM2=ABS(M2)
      MSN=MSS(16)
      MG=ABS(MSS(1))
      MT1=MSS(12)
      MT2=MSS(13)
      MB1=MSS(10)
      MB2=MSS(11)
      MW1=ABS(MSS(27))
      MW2=ABS(MSS(28))
      AMU=ABS(MU)
      MUS=MU**2
      XLAM=DLOG(DBLE(HIGFRZ**2))
C      FTMT=MTMT/V
C      FBMZ=MBMZ/VP
C      FTAMZ=MTAMZ/VP
C          Be careful in using our convention vs Pierce et al.
C          cos(th)>1/sqrt(2):eigenstates same; cs-> -cs
C          cos(th)<1/sqrt(2);flip mass eigenstates; c <-> s interchange
C          Formula remains invariant under these switches
C          Use negative gaugino masses for consistency
C          Now input self-energies for consistency: RSIG*
C
C     here, add in 2 loop QCD correction to mt(DRbar) from Bednyakov
C     et al. Eq. 61.
      ST2LP=CF*(ASMTP/4./PI)**2*(-43.-12*ZETA2+CF*(-59./8.+30*ZETA2-
     ,48*LOG(2.)*ZETA2+12*ZETA3)+
     ,CA*(1093./24.-8*ZETA2+24*LOG(2.)*ZETA2-6*ZETA3))
      MTMT=MT/(1.+5*ASMTP/3./PI+ST2LP)
      FTMT=MTMT/SINB/VEV
      FBMZ=MBMZ/COSB/VEV
      FTAMZ=MTAMZ/COSB/VEV
      RSIGT=SSRSGT(MT**2)
      RSIGB=SSRSGB(MBQ**2)
      RSIGL=SSRSGL(MLQ**2)
C
C Weak Scale Yukawas used by A. Box RGE code
C
      LAMBMZ=DBLE(MBMZ/VEV)
      LAMTAMZ=DBLE(MTAMZ/VEV)
      LAMTMT=DBLE(MTMT/VEV)
C
C     Here, conversion from MSbar to DRbar is done at MZ.
C     Effect of sparticles is included by decoupling
C     beta functions in RGEs
      DAS=DBLE(ASMZ)/2.D0/DPI*(.5)
      ALEM=1./137.036
      DAEM=0.0682-ALEM/2./PI*(-7*LOG(AMW/AMZ))
      ALEMDR=ALEM/(1.-DAEM)
      XWMSB=.23113
C      XW=.2324-1.03E-7*(AMT**2-138.**2)
C     Convert XW to DRbar
      XW=1.-(1.-XWMSB)/(1.-ALEMDR/12./PI)
      A1MZ=5*ALEMDR/3./(1.-XW)
      A2MZ=ALEMDR/XW
C      ALEM=1./128.
C      A1MZ=5*ALEM/3./(1.-XW)
C      A2MZ=ALEM/XW
      G(1)=DSQRT(4*DPI*A1MZ)
      G(2)=DSQRT(4*DPI*A2MZ)
      G(3)=DSQRT(4*DPI*ASMZ/(1.D0-DAS))
      G(4)=DBLE(FTAMZ)
      G(5)=DBLE(FBMZ)
      G(6)=G(6)
      G(25)=DBLE(MU)
      G(26)=DBLE(B)
      G(27)=0.D0
      G(28)=0.D0
      G(29)=0.D0
      G(30)=DBLE(VP)
      G(31)=DBLE(V)
C          Compute gauge mediated threshold functions
      IF (IMODEL.EQ.2) THEN
        XLM=XLAMGM/XMESGM
        THRF=((1.D0+XLM)*(LOG(1.D0+XLM)-2*DDILOG(XLM/(1.D0+XLM))+
     ,        .5*DDILOG(2*XLM/(1.D0+XLM)))+
     ,       (1.D0-XLM)*(LOG(1.D0-XLM)-2*DDILOG(-XLM/(1.D0-XLM))+
     ,        .5*DDILOG(-2*XLM/(1.D0-XLM))))/XLM**2
        THRG=((1.D0+XLM)*LOG(1.D0+XLM)+(1.D0-XLM)*LOG(1.D0-XLM))
     ,/XLM**2
      END IF
C
C          Run back up to mgut with approximate susy spectra
C
      IF (IMODEL.EQ.1) THEN
        IF (XSUGIN(7).EQ.0.) THEN
          MGUT=1.E19
        ELSE
          MGUT=XSUGIN(7)
        END IF
      ELSE IF (IMODEL.EQ.2) THEN
        MGUT=XMESGM
      END IF
      TZ=DLOG(DBLE(MZ)/DBLE(MGUT))
      TGUT=0.D0
      DT=(TGUT-TZ)/DBLE(FLOAT(NSTEP))
      Q=MZ
      DO 250 II=1,NSTEP
        T=TZ+(TGUT-TZ)*FLOAT(II-1)/DBLE(FLOAT(NSTEP))
        Q=SNGL(MGUT*DEXP(T))
        QNEW=SNGL(MGUT*DEXP(T+DT))
        IF (Q.LE.MT.AND.QNEW.GT.MT) G(6)=DBLE(FTMT)
C       Implement sparticle threshold corrections at Q=HIGFRZ
        IF (Q.LE.HIGFRZ.AND.QNEW.GT.HIGFRZ) THEN
          G(6)=G(6)/(1.D0-DBLE(RSIGT))
          G(5)=G(5)/(1.D0-DBLE(RSIGB))
          G(4)=G(4)/(1.D0-DBLE(RSIGL))
          IF (INUHM.EQ.1) THEN
            G(13)=DBLE(MHDSQ)
            G(14)=DBLE(MHUSQ)
          END IF
        END IF
        IF (Q.LE.XNRIN(2).AND.QNEW.GT.XNRIN(2)) THEN
          G(27)=DBLE(FNMZ)
          G(28)=DBLE(G0(28))
          G(29)=DBLE(G0(29))
        END IF
        CALL DRKSTP(31,DT,T,G,SURG26,W2)
        A1I=SNGL(4*DPI/G(1)**2)
        A2I=SNGL(4*DPI/G(2)**2)
        A3I=SNGL(4*DPI/G(3)**2)
C       TEST YUKAWA DIVERGENCE
        IF (G(4).GT.5.D0.OR.G(5).GT.5.D0.OR.
     $G(6).GT.5.D0.OR.G(27).GT.5.D0) THEN
          NOGOOD=4
          GO TO 100
        END IF
        IF (A1I.LT.A2I.AND.XSUGIN(7).EQ.0.) GO TO 30
250   CONTINUE
      IF (IMODEL.EQ.1.AND.XSUGIN(7).EQ.0.) THEN
        WRITE(LOUT,*) 'SUGRGE ERROR: NO UNIFICATION FOUND'
        NOGOOD=1
        GO TO 100
      END IF
30    IF (XSUGIN(7).EQ.0.) THEN
        MGUT=QNEW
      ELSE
        MGUT=XSUGIN(7)
      END IF
      AGUT=SNGL((G(1)**2/4.D0/DPI+G(2)**2/4.D0/DPI)/2.D0)
      GGUT=SQRT(4*PI*AGUT)
      AGUTI=1./AGUT
      FTAGUT=SNGL(G(4))
      FBGUT=SNGL(G(5))
      FTGUT=SNGL(G(6))
      IF (INUHM.EQ.1) THEN
        MHDSMG=SNGL(G(13))
        MHUSMG=SNGL(G(14))
      END IF
      MUMG=SNGL(G(25))
      BMG=SNGL(G(26))
      IF (XNRIN(2).LT.1.E19.AND.XNRIN(1).EQ.0.) THEN
C     IMPOSE FN-FT UNIFICATION
        FNGUT=SNGL(G(6))
      ELSE
        FNGUT=SNGL(G(27))
      END IF
      G3GUT=SNGL(G(3))
      MGUTSS=MGUT
      AGUTSS=AGUT
      GGUTSS=GGUT
C
C          Set GUT boundary condition
C
      DO 260 I=1,3
        IF (IMODEL.EQ.1) THEN
          G(I)=G(I)
          G(I+6)=DBLE(MHF)
          G(I+9)=DBLE(A0)
        ELSE IF (IMODEL.EQ.2) THEN
          G(I)=G(I)
          G(I+6)=XGMIN(11+I)*XGMIN(8)*THRG*(G(I)/4.D0/DPI)**2*XLAMGM
          G(I+9)=0.D0
        END IF
      IF (XNRIN(2).LT.1.E19) THEN
        G(27)=DBLE(FNGUT)
        G(28)=DBLE(XNRIN(4))**2
        G(29)=DBLE(XNRIN(3))
      ELSE
        G(27)=0.D0
        G(28)=0.D0
        G(29)=0.D0
      END IF
260   CONTINUE
C     OVERWRITE ALFA_3 UNIFICATION TO GET ALFA_3(MZ) RIGHT
      IF (IMODEL.EQ.1.AND.IAL3UN.NE.0) G(3)=DBLE(GGUT)
      IF (IMODEL.EQ.1) THEN
        DO 270 I=13,24
          G(I)=DBLE(M0)**2
270     CONTINUE
      IF (INUHM.EQ.1) THEN
        G(13)=DBLE(MHDSMG)
        G(14)=DBLE(MHUSMG)
      END IF
C          Set possible non-universal GUT scale boundary conditions
      DO 280 I=1,6
        IF (XNUSUG(I).LT.1.E19) THEN
          G(I+6)=DBLE(XNUSUG(I))
        END IF
280   CONTINUE
      DO 281 I=7,18
        IF (XNUSUG(I).LT.1.E19) THEN
          G(I+6)=SIGN(1.,XNUSUG(I))*DBLE(XNUSUG(I))**2
        END IF
281   CONTINUE
      ELSE IF (IMODEL.EQ.2) THEN
       XC=2*THRF*XLAMGM**2
       DY=DSQRT(3.D0/5.D0)*G(1)*XGMIN(11)
       G(13)=XC*(.75*XGMIN(13)*(G(2)/4.D0/DPI)**4+.6D0*.25*
     , XGMIN(12)*(G(1)/4.D0/DPI)**4)+DBLE(XGMIN(9))-DY
       G(14)=XC*(.75*XGMIN(13)*(G(2)/4.D0/DPI)**4+.6D0*.25*
     , XGMIN(12)*(G(1)/4.D0/DPI)**4)+DBLE(XGMIN(10))+DY
       G(15)=XC*(.6*XGMIN(12)*(G(1)/4.D0/DPI)**4)+2*DY
       G(16)=XC*(.75*XGMIN(13)*(G(2)/4.D0/DPI)**4+.6D0*.25*
     , XGMIN(12)*(G(1)/4.D0/DPI)**4)-DY
       G(17)=XC*(4*XGMIN(14)*(G(3)/4.D0/DPI)**4/3.D0+.6D0*XGMIN(12)*
     , (G(1)/4.D0/DPI)**4/9.D0)+2*DY/3.D0
       G(18)=XC*(4*XGMIN(14)*(G(3)/4.D0/DPI)**4/3.D0+
     , .6D0*4*XGMIN(12)*(G(1)/4.D0/DPI)**4/9.D0)-4*DY/3.D0
       G(19)=XC*(4*XGMIN(14)*(G(3)/4.D0/DPI)**4/3.D0+.75*XGMIN(13)*
     ,(G(2)/4.D0/DPI)**4+.6*XGMIN(12)*(G(1)/4.D0/DPI)**4/36.D0)
     ,+DY/3.D0
       G(20)=G(15)
       G(21)=G(16)
       G(22)=G(17)
       G(23)=G(18)
       G(24)=G(19)
      ELSE IF (IMODEL.EQ.7.OR.IMODEL.EQ.9.OR.IMODEL.EQ.10) THEN
       G(1)=G(1)
       G(2)=G(2)
       G(3)=G(3)
       BLHAT=G(4)*(-9*G(1)**2/5.D0-3*G(2)**2+3*G(5)**2+4*G(4)**2)
       BBHAT=G(5)*(-7*G(1)**2/15.D0-3*G(2)**2-16*G(3)**2/3.D0+
     ,             G(6)**2+6*G(5)**2+G(4)**2)
       BTHAT=G(6)*(-13*G(1)**2/15.D0-3*G(2)**2-16*G(3)**2/3.D0+
     ,             6*G(6)**2+G(5)**2)
       G(7)=33*MHF*G(1)**2/5.D0/16.D0/DPI**2
       IF (IMODEL.EQ.10) THEN
         G(7)=G(7)+XAMIN(11)*MHF
       END IF
       G(8)=MHF*G(2)**2/16.D0/DPI**2
       G(9)=-3*MHF*G(3)**2/16.D0/DPI**2
       G(10)=-BLHAT*MHF/G(4)/16.D0/DPI**2
       G(11)=-BBHAT*MHF/G(5)/16.D0/DPI**2
       G(12)=-BTHAT*MHF/G(6)/16.D0/DPI**2
       G(13)=(-99*G(1)**4/50.D0-3*G(2)**4/2.D0+3*G(5)*BBHAT+
     ,G(4)*BLHAT)*MHF**2/(16*DPI**2)**2+XAMIN(6)*DBLE(M0)**2
       G(14)=(-99*G(1)**4/50.D0-3*G(2)**4/2.D0+3*G(6)*BTHAT)*
     ,        MHF**2/(16*DPI**2)**2+XAMIN(7)*DBLE(M0)**2
       G(15)=(-198*G(1)**4/25.D0)*MHF**2/(16*DPI**2)**2+
     ,XAMIN(5)*DBLE(M0)**2
       G(16)=(-99*G(1)**4/50.D0-3*G(2)**4/2.D0)*MHF**2/(16*DPI**2)**2
     ,+XAMIN(4)*DBLE(M0)**2
       G(17)=(-22*G(1)**4/25.D0+8*G(3)**4)*MHF**2/(16*DPI**2)**2+
     ,XAMIN(2)*DBLE(M0)**2
       G(18)=(-88*G(1)**4/25.D0+8*G(3)**4)*MHF**2/(16*DPI**2)**2+
     ,XAMIN(3)*DBLE(M0)**2
       G(19)=(-11*G(1)**4/50.D0-3*G(2)**4/2.D0+8*G(3)**4)*
     ,        MHF**2/(16*DPI**2)**2+XAMIN(1)*DBLE(M0)**2
       G(20)=(-198*G(1)**4/25.D0+2*G(4)*BLHAT)*MHF**2/(16*DPI**2)**2
     ,+XAMIN(5)*DBLE(M0)**2
       G(21)=(-99*G(1)**4/50.D0-3*G(2)**4/2.D0+G(4)*BLHAT)*
     ,        MHF**2/(16*DPI**2)**2+XAMIN(4)*DBLE(M0)**2
       G(22)=(-22*G(1)**4/25.D0+8*G(3)**4+2*G(5)*BBHAT)*
     , MHF**2/(16*DPI**2)**2+XAMIN(2)*DBLE(M0)**2
       G(23)=(-88*G(1)**4/25.D0+8*G(3)**4+2*G(6)*BTHAT)*
     , MHF**2/(16*DPI**2)**2+XAMIN(3)*DBLE(M0)**2
       G(24)=(-11*G(1)**4/50.D0-3*G(2)**4/2.D0+8*G(3)**4+G(5)*BBHAT+
     ,        G(6)*BTHAT)*MHF**2/(16*DPI**2)**2+XAMIN(1)*DBLE(M0)**2
      END IF
      IF (IMODEL.EQ.9) THEN
        CALL MMAMSB(M0,MHF,G)
      END IF
      DO 285 I=1,31
        IG(I)=0
285   CONTINUE
C          Check for tachyonic sleptons at GUT scale
      IF (G(15).LT.0.D0.OR.G(16).LT.0.D0) THEN
        ITACHY=2
      ELSE
        ITACHY=0
      END IF
C
C          Run back down to weak scale
C
      TZ=DLOG(DBLE(MZ)/DBLE(MGUT))
      TGUT=0.D0
      DT=(TZ-TGUT)/DBLE(FLOAT(NSTEP))

      DO II=1,10000
         TQV(II)=0d0
         TRHM1V(II)=0d0
         TRHM2V(II)=0d0
         TRHM3V(II)=0d0
         TRHG1V(II)=0d0
         TRHG2V(II)=0d0
         TRHG3V(II)=0d0
      ENDDO


      DO 290 II=1,NSTEP+2
        T=TGUT+(TZ-TGUT)*FLOAT(II-1)/DBLE(FLOAT(NSTEP))
        Q=SNGL(MGUT*DEXP(T))
        CALL DRKSTP(31,DT,T,G,SURG26,W2)
C       Here, DRKSTP advances T by DT
        QNEW=SNGL(MGUT*DEXP(T))
C       TEST YUKAWA DIVERGENCE
        IF (G(4).GT.5.D0.OR.G(5).GT.5.D0.OR.
     $    G(6).GT.5.D0.OR.G(27).GT.5.D0) THEN
          NOGOOD=4
          GO TO 100
        END IF
        CALL SUGFRZ(QNEW,G,G0,IG)
        IF (Q.GE.AMNRMJ.AND.QNEW.LT.AMNRMJ.AND.XNRIN(1).EQ.0.) THEN
          FNMZ=SNGL(G(27))
        END IF
        IF (Q.GT.HIGFRZ.AND.QNEW.LE.HIGFRZ) THEN
          G(6)=G(6)*(1.D0-DBLE(RSIGT))
          G(5)=G(5)*(1.D0-DBLE(RSIGB))
          G(4)=G(4)*(1.D0-DBLE(RSIGL))
        END IF
        IF (QNEW.LT.AMNRMJ) THEN
          G(27)=0.D0
          G(28)=0.D0
          G(29)=0.D0
        END IF
        IF (NOGOOD.NE.0) GO TO 100
        IF (QNEW.LT.MZ) GO TO 40

C...@@@ SAVE RUNNING PARAMETERS:    
      IF(II.GT.100000) THEN
         WRITE(*,*) 'ERROR: TOO MANY POINTS IN RGE'
         STOP
      ELSE
         TQV(II)=DBLE(Q)
         TRHM1V(II)=DBLE(G(7))
         TRHM2V(II)=DBLE(G(8))
         TRHM3V(II)=DBLE(G(9))
         TRHG1V(II)=DBLE(G(1))
         TRHG2V(II)=DBLE(G(2))
         TRHG3V(II)=DBLE(G(3))
      ENDIF



290   CONTINUE
40    CONTINUE
C
C          Electroweak breaking constraints; tree level
C
      VUQ=G0(31)
      VDQ=G0(30)
      TANBQ=VUQ/VDQ
      SIN2BQ=SIN(2*ATAN(TANBQ))
      MZQ=SQRT((G0(2)**2+.6*G0(1)**2)*(VUQ**2+VDQ**2)/2.)
      BADMU=.FALSE.
      IF (INUHM.NE.1) THEN
      MUS=(G0(13)-G0(14)*TANBQ**2)/(TANBQ**2-1.)-MZQ**2/2.
C          Compute loop corrections using masses from last iteration
      CALL SUGEFF(G0,SIG1,SIG2,SIG3)
      MH1S=G0(13)+SIG1
      MH2S=G0(14)+SIG2
      MUS=(MH1S-MH2S*TANBQ**2)/(TANBQ**2-1.)-MZQ**2/2.
C          If MUS<0, set it to MZ**2 and continue
      IF (MUS.LT.0.) THEN
        MUS=AMZ**2
      END IF
      MU=SQRT(MUS)*SIGN(1.,SGNMU)
      B=(G0(13)+G0(14)+2*MUS)*SIN2BQ/MU/2.+SIG3/MU
      CALL SUGMAS(G0,0,IMODEL,SIGA)
      IF (NOGOOD.NE.0) GO TO 100
C
C           Electroweak breaking constraints; loop level
C
      CALL SUGEFF(G0,SIG1,SIG2,SIG3)
      MH1S=G0(13)+SIG1
      MH2S=G0(14)+SIG2
      MUS=(MH1S-MH2S*TANBQ**2)/(TANBQ**2-1.)-MZQ**2/2.
      IF (MUS.LT.0.) THEN
C        NOGOOD=2
C        GO TO 100
         MUS=MZ**2
      END IF
      MU=SQRT(MUS)*SIGN(1.,SGNMU)
      B=(MH1S+MH2S+2*MUS)*SIN2BQ/MU/2.+SIG3/MU
C
C     Once more, with feeling!
C
      CALL SUGEFF(G0,SIG1,SIG2,SIG3)
      MH1S=G0(13)+SIG1
      MH2S=G0(14)+SIG2
      MUS=(MH1S-MH2S*TANBQ**2)/(TANBQ**2-1.)-MZQ**2/2.
      IF (MUS.LT.0.) THEN
C        NOGOOD=2
C        GO TO 100
         BADMU=.TRUE.
         MUS=MZ**2
      END IF
      MU=SQRT(MUS)*SIGN(1.,SGNMU)
      B=(MH1S+MH2S+2*MUS)*SIN2BQ/MU/2.+SIG3/MU
      CALL SUGMAS(G0,1,IMODEL,SIGA)
      ELSE
        MUS=MU**2
        B=AMHA**2/MU/(COTB+TANB)+SIG3/MU
        CALL SUGMAS(G0,1,IMODEL,SIGA)
        CALL SUGEFF(G0,SIG1,SIG2,SIG3)
        MHDSQ=(AMHA**2+MZQ**2)*TANB**2/(TANB**2+1.)-SIG1-MUS-MZQ**2/2.
      MHUSQ=(AMHA**2+MZQ**2)/(TANB**2+1.)-MUS-MZQ**2/2.-SIG2
      END IF
C
C  Save radiative corrections to Yukawas for b->s gamma computation
C
      RTISA=RSIGT
      RBISA=RSIGB
      RLISA=RSIGL
C
C Initial value of MU passed to A. Box RGE code
C
      RGEMU=-MU !THE -VE SIGN FIXES CONVENTIONS
C

100   RETURN
      END














































CDECK  ID>, SSPRT.  
      SUBROUTINE PLOTBF(MOM,DAUG,RAT,BF)
C-----------------------------------------------------------------------
C
C     Returns BF for ID1->ID2
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
      COMMON/SSLUN/LOUT,LHEOUT
      INTEGER LOUT,LHEOUT
      SAVE /SSLUN/
C          MXSS         =  maximum number of modes
C          NSSMOD       = number of modes
C          ISSMOD       = initial particle
C          JSSMOD       = final particles
C          GSSMOD       = width
C          BSSMOD       = branching ratio
C          MSSMOD       = decay matrix element pointer
C          LSSMOD       = logical flag used internally by SSME3
      INTEGER MXSS
      PARAMETER (MXSS=1000)
      COMMON/SSMODE/NSSMOD,ISSMOD(MXSS),JSSMOD(5,MXSS),GSSMOD(MXSS)
     $,BSSMOD(MXSS),MSSMOD(MXSS),LSSMOD
      INTEGER NSSMOD,ISSMOD,JSSMOD,MSSMOD
      REAL GSSMOD,BSSMOD
      LOGICAL LSSMOD
      SAVE /SSMODE/
C

      INTEGER ID1,ID2,I,K,NOUT
      REAL BF,RAT
      CHARACTER*5 MOM,DAUG(4)
      CHARACTER*5 SSID,LBLIN,LBLOUT(4)
C

      NOUT=0
      DO 100 I=1,NSSMOD
        IF(SSID(ISSMOD(I)).NE.MOM) GOTO 100
        NOUT=NOUT+1
        DO 110 K=1,4
          LBLOUT(K)=SSID(JSSMOD(K,I))
          IF(LBLOUT(K).NE.DAUG(K)) GOTO 100
110     CONTINUE
        RAT=GSSMOD(I)
        BF=BSSMOD(I)
        GOTO 111
100   CONTINUE
      BF=0.
      RAT=0.      ! Set to zero if mode not found
111   RETURN
      END

CDECK  ID>, SSID.   
      CHARACTER*5 FUNCTION SSID(ID)
C-----------------------------------------------------------------------
C
C     Return character name for ID, assuming the default IDENT codes
C     are used in /SSTYPE/.
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
      COMMON/SSLUN/LOUT,LHEOUT
      INTEGER LOUT,LHEOUT
      SAVE /SSLUN/
      CHARACTER*5 LABEL(-120:120)
      SAVE LABEL
      INTEGER ID,J
C
      DATA LABEL(0)/'     '/
C
      DATA (LABEL(J),J=1,10)
     $/'UP   ','DN   ','ST   ','CH   ','BT   ','TP   '
     $,'ERROR','ERROR','GL   ','GM   '/
      DATA (LABEL(J),J=-1,-10,-1)
     $/'UB   ','DB   ','SB   ','CB   ','BB   ','TB   '
     $,'ERROR','ERROR','ERROR','ERROR'/
C
      DATA (LABEL(J),J=11,20)
     $/'NUE  ','E-   ','NUM  ','MU-  ','NUT  ','TAU- '
     $,'ERROR','ERROR','ERROR','ERROR'/
      DATA (LABEL(J),J=-11,-20,-1)
     $/'ANUE ','E+   ','ANUM ','MU+  ','ANUT ','TAU+ '
     $,'ERROR','ERROR','ERROR','ERROR'/
C
      DATA (LABEL(J),J=21,30)
     $/'UPL  ','DNL  ','STL  ','CHL  ','BT1  ','TP1  '
     $,'ERROR','ERROR','GLSS ','Z1SS '/
      DATA (LABEL(J),J=-21,-30,-1)
     $/'UBL  ','DBL  ','SBL  ','CBL  ','BB1  ','TB1  '
     $,'ERROR','ERROR','ERROR','ERROR'/
C
      DATA (LABEL(J),J=31,40)
     $/'NUEL ','EL-  ','NUML ','MUL- ','NUTL ','TAU1-'
     $,'ERROR','ERROR','W1SS+','Z2SS '/
      DATA (LABEL(J),J=-31,-40,-1)
     $/'ANUEL','EL+  ','ANUML','MUL+ ','ANUTL','TAU1+'
     $,'ERROR','ERROR','W1SS-','ERROR'/
C
      DATA (LABEL(J),J=41,50)
     $/'UPR  ','DNR  ','STR  ','CHR  ','BT2  ','TP2  '
     $,'ERROR','ERROR','W2SS+','Z3SS '/
      DATA (LABEL(J),J=-41,-50,-1)
     $/'UBR  ','DBR  ','SBR  ','CBR  ','BB2  ','TB2  '
     $,'ERROR','ERROR','W2SS-','ERROR'/
C
      DATA (LABEL(J),J=51,60)
     $/'NUER ','ER-  ','NUMR ','MUR- ','NUTR ','TAU2-'
     $,'ERROR','ERROR','ERROR','Z4SS '/
      DATA (LABEL(J),J=-51,-60,-1)
     $/'ANUEL','ER+  ','ANUMR','MUR+ ','ANUTR','TAU2+'
     $,'ERROR','ERROR','ERROR','ERROR'/
C
      DATA (LABEL(J),J=82,86)
     $/'HL0  ','HH0  ','HA0  ','ERROR','H+   '/
      DATA LABEL(-86)/'H-   '/
C
      DATA LABEL(80)/'W+   '/,LABEL(-80)/'W-   '/,LABEL(90)/'Z0   '/
      DATA LABEL(91)/'GVSS '/
      DATA LABEL(110)/'PI0  '/
      DATA LABEL(120)/'PI+  '/,LABEL(-120)/'PI-  '/
C

      IF(IABS(ID).GT.120) THEN
        WRITE(LOUT,*) 'SSID: ID = ',ID
        STOP99
      ENDIF
      SSID=LABEL(ID)
      RETURN
      END

C...#######################################

      SUBROUTINE FINDLSP(ID)
C...Check all masses and return the ID of the LSP

      IMPLICIT NONE


C          Masses:
C     MSS( 1) = glss     MSS( 2) = upl      MSS( 3) = upr
C     MSS( 4) = dnl      MSS( 5) = dnr      MSS( 6) = stl
C     MSS( 7) = str      MSS( 8) = chl      MSS( 9) = chr
C     MSS(10) = b1       MSS(11) = b2       MSS(12) = t1
C     MSS(13) = t2       MSS(14) = nuel     MSS(15) = numl
C     MSS(16) = nutl     MSS(17) = el-      MSS(18) = er-
C     MSS(19) = mul-     MSS(20) = mur-     MSS(21) = tau1
C     MSS(22) = tau2     MSS(23) = z1ss     MSS(24) = z2ss
C     MSS(25) = z3ss     MSS(26) = z4ss     MSS(27) = w1ss
C     MSS(28) = w2ss     MSS(29) = hl0      MSS(30) = hh0
C     MSS(31) = ha0      MSS(32) = h+
C          Unification:
C     MGUTSS  = M_GUT    GGUTSS  = g_GUT    AGUTSS  = alpha_GUT
      COMMON /SUGMG/ MSS(32),GSS(31),MGUTSS,GGUTSS,AGUTSS,FTGUT,
     $FBGUT,FTAGUT,FNGUT
      REAL MSS,GSS,MGUTSS,GGUTSS,AGUTSS,FTGUT,FBGUT,FTAGUT,FNGUT

      INTEGER ID,I
      REAL MLSP

      MLSP=ABS(MSS(1))
      DO I=1,28
         IF(ABS(MSS(I)).LT.MLSP) THEN
            MLSP=ABS(MSS(I))
            ID=I
         ENDIF
      ENDDO
       
      END

      DOUBLE PRECISION FUNCTION PSLAMB(X,Y,Z)
      
      IMPLICIT NONE
      
      DOUBLE PRECISION X,Y,Z
      
      PSLAMB = X**2 + Y**2 + Z**2 - 2d0*X*Y - 2d0*X*Z - 2d0*Z*Y
      
      END FUNCTION
      
      
      DOUBLE PRECISION FUNCTION COT(X)
      
      IMPLICIT NONE
      
      DOUBLE PRECISION X
      
      COT = 1d0/DTAN(X)
      
      END FUNCTION
