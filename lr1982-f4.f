C      This code is taken from "A Guide to the Use of the ITS Irregular
C      Terrain Model in the Area Prediction Mode," by G.A. Hufford,
C      A.G. Longley, and W.A. Kissick.  It is ITM v1.2.1 from 1982.
C      provided for historical interest and superceded by newer
C      versions.  Please do not use this for your propagation
C      predictions.  This does not include the 1985 corrections by
C      Hufford.
C
C      You will likely need a FORTRAN-IV compiler like ifort.  While it
C      compiles under f77 with warnings, results include some NaNs.
C      After compiling, run without command line arguments and output
C      will be be printed to stdout.  For example,
C
C      ifort -o lr1982 lr1982-f4.f
C      lr1982 < inputDeck | more
C
C      File inputDeck follows, where of course leading 'C' is removed,
C      but not leading spaces:
C=== inputDeck follows ===
C 1    
CTEST PROBLEM 1, QKAREA
C 2        0.        100.      10.       400.      25.
C 32       70.
C 4        10.       50.       90.       95.
C870110    400.      10.       1.        200.
C 1
CTEST PROBLEM 2, QKAREA
C861       25.
C99
C 1
CTEST PROBLEM 4, QKAREA
C 31       10.
C 57       30.       350.      300.      25.       0.02
C 60000    1200.     5.        1000.
C88
C
C=== inputDeck ends after required blank line ===
C
C      Mike Markowski, mike.ab3ap@gmail.com
C      Dec 2022

      SUBROUTINE LRPROP(D)
C        COMPUTES AREF, THE REFERENCE VALUE OF RADIO ATTENUATION
C          VERSION 1.2.1 (AUG 71/MAR 77/APR 79)
C             OF THE LONGLEY-RICE (1968) MODEL
C          PRINCIPAL CHANGES-
C             1.1. A SIMPLIFICATION OF THE LINE-OF-SIGHT AND SCATTER
C                  ROUTINES
C             1.2. A CHANGE IN THE LINE-OF-SIGHT ROUTINE AND IN THE
C                  SUBSEQUENT CALCULATIONS. RESULTS ARE IMPROVED WHEN
C                  ONE OR BOTH ANTENNAS ARE HIGH.
C          VALID ONLY FOR...
C               FREQUENCIES BETWEEN 20 MHZ AND 20 GHZ
C               ANTENNA HEIGHTS BETWEEN 0.5 M AND 3000 M
C               ELEVATION ANGLES LESS THAN 200 MRAD
C
      COMMON/PROP/KWX,AREF,MDP,DIST,HG(2),WN,DH,ENS,GME,ZGND,
     X   HE(2),DL(2),THE(2)
         COMPLEX ZGND
      COMMON/PROPA/DLSA,DX,AEL,AK1,AK2,AED,EMD,AES,EMS,DLS(2),DLA,THA
C
      COMMON/SAVE/SAVA(6),WLOS,WSCAT,DMIN,XAE,SAVB(40)
C
      LOGICAL WLOS,WSCAT
C
      DATA THIRD/0.3333333/
C
       IF(MDP) 10,32,10
C
  10   CONTINUE
       DO 11 J=1,2
  11   DLS(J)=SQRT(2.*HE(J)/GME)
      DLSA=DLS(1)+DLS(2)
      DLA=DL(1)+DL(2)
      THA=AMAX1(THE(1)+THE(2),-DLA*GME)
      WLOS=.FALSE.
      WSCAT=.FALSE.
C
C       CHECK PARAMETER RANGES
         IF(ENS .LT. 250. .OR. ENS .GT. 400.) GO TO 154
         IF(GME .LT. 75E-9 .OR. GME .GT. 250E-9) GO TO 154
         IF(REAL(ZGND) .LE. ABS(AIMAG(ZGND))) GO TO 154
         DO 121 J=1,2
           IF(ABS(THE(J)) .GT. 200E-3) GO TO 153
           IF(DL(J) .LT. 0.1*DLS(J) .OR. DL(J) .GT. 3.*DLS(J))
     X       GO TO 153
  121    CONTINUE
         IF(WN .LT. 0.838 .OR. WN .GT. 210.) GO TO 151
         DO 122 J=1,2
           IF(HG(J) .LT. 1. .OR. HG(J) .GT. 1000.) GO TO 151
  122    CONTINUE
       GO TO 158
  153     KWX=MAX0(KWX,3)
  151     KWX=MAX0(KWX,1)
          IF(WN .LT. 0.419 .OR. WN .GT. 420.) GO TO 154
          DO 132 J=1,2
            IF(HG(J) .LT. 0.5 .OR. HG(J) .GT. 3000.) GO TO 154
  132     CONTINUE
        GO TO 158
  154     KWX=4
  158  CONTINUE
       DMIN=ABS(HE(1)-HE(2))/200E-3
C
C        COEFFICIENTS FOR THE DIFFRACTION RANGE
C
       Q=ADIFF(0.)
       XAE=(WN*GME**2)**(-THIRD)
       D3=AMAX1(DLSA,1.3787*XAE+DLA)
       D4=D3+2.7574*XAE
       A3=ADIFF(D3)
       A4=ADIFF(D4)
       EMD=(A4-A3)/(D4-D3)
       AED=A3-EMD*D3
C
       IF(MDP) 33,32,31
  31      MDP=0
  32      DIST=D
        IF(DIST .LE. 0.) GO TO 38
  33      CONTINUE
          IF(DIST .GT. 1000E3) KWX=MAX0(KWX,1)
          IF(DIST .LT. DMIN) KWX=MAX0(KWX,3)
          IF(DIST .LT. 1E3 .OR. DIST .GT. 2000E3) KWX=4
  38    CONTINUE
C
        IF(DIST .GE. DLSA) GO TO 50
C
        IF(WLOS) GO TO 48
C
C         COEFFICIENTS FOR THE LINE-OF-SIGHT RANGE
C
      Q=ALOS(0.)
      D2=DLSA
      A2=AED+D2*EMD
      D0=1.908*WN*HE(1)*HE(2)
      IF(AED .LT. 0.) GO TO 41
         D0=AMIN1(D0,0.5*DLA)
         D1=D0+0.25*(DLA-D0)
       GO TO 42
  41     D1=AMAX1(-AED/EMD,0.25*DLA)
  42  A1=ALOS(D1)
      IF(D0 .GE. D1) GO TO 43
         A0=ALOS(D0)
         Q=ALOG(D2/D0)
         AK2=AMAX1(0.,((D2-D0)*(A1-A0)-(D1-D0)*(A2-A0))/
     X      ((D2-D0)*ALOG(D1/D0)-(D1-D0)*Q))
       IF(AK2 .GT. 0.) GO TO 44
       IF(AED .GE. 0.) GO TO 44
  43      AK2=0.
          AK1=(A2-A1)/(D2-D1)
          IF(AK1 .GT. 0.) GO TO 46
          GO TO 45
  44         AK1=(A2-A0-AK2*Q)/(D2-D0)
             IF(AK1 .GE. 0.) GO TO 46
             AK1=0.
             AK2=DIM(A2,A0)/Q
           IF(AK2 .GT. 0.) GO TO 46
  45         AK1=EMD
  46      AEL=A2-AK1*D2-AK2*ALOG(D2)
          WLOS=.TRUE.
  48      IF(DIST .LE. 0.) GO TO 50
          AREF=AEL+AK1*DIST+AK2*ALOG(DIST)
          GO TO 60
C
  50      IF(WSCAT) GO TO 58
C
C           COEFFICIENTS FOR THE SCATTER RANGE
C
          Q=ASCAT(0.)
          D5=DLA+200E3
          D6=D5+200E3
          A6=ASCAT(D6)
          A5=ASCAT(D5)
          IF(A5 .LT. 1000.) GO TO 51
             EMS=EMD
             AES=AED
             DX=10E6
           GO TO 52
  51         EMS=(A6-A5)/200E3
             DX=AMAX1(DLSA,DLA+0.3*XAE*ALOG(47.7*WN),
     X          (A5-AED-EMS*D5)/(EMD-EMS))
             AES=(EMD-EMS)*DX+AED
  52      WSCAT=.TRUE.
C
  58      IF(DIST .GT. DX) GO TO 59
             AREF=AED+EMD*DIST
           GO TO 60
  59         AREF=AES+EMS*DIST
C
  60      AREF=DIM(AREF,0.)
          RETURN
         END

      FUNCTION ADIFF(D)
C        THE *DIFFRACTION ATTENUATION* AT DISTANCE D
C        A CONVEX COMBINATION OF SMOOTH EARTH DIFFRACTION AND
C           DOUBLE KNIFE-EDGE DIFFRACTION
C        A CALL WITH D=0 SETS UP INITIAL CONSTANTS
C
      COMMON/PROP/KWX,AREF,MDP,DIST,HG(2),WN,DH,ENS,GME,ZGND,
     X   HE(2),DL(2),THE(2)
      COMPLEX ZGND
      COMMON/PROPA/DLSA,DX,AEL,AK1,AK2,AED,EMD,AES,EMS,DLS(2),DLA,THA
C
      COMMON/SAVE/WD1,XD1,AFO,QK,AHT,XHT,SAVE(44)
C
      DATA THIRD/0.3333333/
C
       IF(D .GT. 0.) GO TO 10
C
       Q=HG( 1 )*HG(2)
       QK=HE(1)*HE(2)-Q
       IF(MDP .LT. 0) Q=Q+10.
       WD1=SQRT(1.+QK/Q)
       XD1=DLA+THA/GME
       Q=(1.-0.5*EXP(-DLSA/50E3))*DH
       Q=0.78*Q*EXP(-(Q/16.)**0.25)
       AFO=AMIN1(15.,2.171*ALOG(1.+4.77E-4*HG(1)*HG(2)*WN*Q))
       QK=1./CABS(ZGND)
       AHT=20.
       XHT=0.
       DO 1 J=1,2
         A=0.5*DL(J)**2/HE(J)
         WA=(A*WN)**THIRD
         PK=QK/WA
         Q=(1.607-PK)*151.0*WA*DL(J)/A
         XHT=XHT+Q
         AHT=AHT+FHT(Q,PK)
   1   CONTINUE
       ADIFF=0.
       GO TO 80
C
  10   CONTINUE
       TH=THA+D*GME
       DS=D-DLA
       Q=0.0795775*WN*DS*TH**2
       ADIFF=AKNFE(Q*DL(1)/(DS+DL(1)))+AKNFE(Q*DL(2)/(DS+DL(2)))
       A=DS/TH
       WA=(A*WN)**THIRD
       PK=QK/WA
       Q=(1.607-PK)*151.0*WA*TH+XHT
       AR=0.05751*Q-4.343*ALOG(Q)-AHT
       Q=(WD1+XD1/D)*AMIN1(((1.-0.8*EXP(-D/50E3))*DH*WN),6283.2)
       WD=25.1/(25.1+SQRT(Q))
       ADIFF=(AR-ADIFF)*WD+ADIFF+AFO
  80   RETURN
      END

      FUNCTION ALOS(D)
C        THE *LINE-OF-SIGHT ATTENUATION* AT DISTANCE D
C        A CONVEX COMBINATION OF PLANE EARTH FIELDS AND
C           DIFFRACTED FIELDS
C        A CALL WITH D=0 SETS UP INITIAL CONSTANTS
C
      COMMON/PROP/KWX,AREF,MDP,DIST,HG(2),WN,DH,ENS,GME,ZGND,
     X   HE(2),DL(2),THE(2)
      COMPLEX ZGND
      COMMON/PROPA/DLSA,DX,AEL,AK1,AK2,AED,EMD,AES,EMS,DLS(2),DLA,THA
C
      COMMON/SAVE/WLS,SAVE(49)
C
      COMPLEX R
C
      ABQ(R)=REAL(R)**2+AIMAG(R)**2
C
       IF(D .GT. 0.) GO TO 10
C
       WLS=0.021/(0.021+WN*DH/AMAX1(10E3,DLSA))
       ALOS=0.
       GO TO 80
C
   10  CONTINUE
       Q=(1.-0.8*EXP(-D/50E3))*DH
       S=0.78*Q*EXP(-(Q/16.)**0.25)
       Q=HE(1)+HE(2)
       SPS=Q/SQRT(D**2+Q**2)
       R=(SPS-ZGND)/(SPS+ZGND)*EXP(-WN*S*SPS)
       Q=ABQ(R)
       IF(Q .LT. 0.25 .OR. Q .LT. SPS) R=R*SQRT(SPS/Q)
       ALOS=EMD*D+AED
       Q=WN*HE(1)*HE(2)*2./D
       ALOS=(-4.343*ALOG(ABQ(CMPLX(COS(Q),-SIN(Q))+R))-ALOS)*WLS+ALOS
C
   80  RETURN
      END

      FUNCTION ASCAT(D)
C        THE *SCATTER ATTENUATION* AT DISTANCE D
C        USES AN APPROXIMATION TO THE METH0DS OF NBS TN101 WITH
C           CHECKS FOR INADMISSABLE SITUATIONS
C        FOR PROPER OPERATION, THE LARGER DISTANCE (D=D6)
C           MUST BE THE FIRST CALLED
C        A CALL WITH D=0 SETS UP INITIAL CONSTANTS
C
      COMMON/PROP/KWX,AREF,MDP,DIST,HG(2),WN,DH,ENS,GME,ZGND,
     X   HE(2),DL(2),THE(2)
      COMPLEX ZGND
      COMMON/PROPA/DLSA,DX,AEL,AK1,AK2,AED,EMD,AES,EMS,DLS(2),DLA,THA
C
      COMMON/SAVE/AD,RR,ETQ,H0S,SAVE(46)
C
       IF(D .GT. 0.) GO TO 10
C
       AD=DL(1)-DL(2)
       RR=HE(2)/HE(1)
       IF(AD) 1,2,2
  1       AD=-AD
          RR=1./RR
  2    ETQ=(5.67E-6*ENS-2.32E-3)*ENS+0.031
       H0S=-15.
       ASCAT=0.
       GO TO 80
C
   10  CONTINUE
       IF(H0S .GT. 15.) GO TO 12
       TH=THE(1)+THE(2)+D*GME
       R2=2.*WN*TH
       R1=R2*HE(1)
       R2=R2*HE(2)
       IF(R1 .GT. 0.2 .OR. R2 .GT. 0.2) GO TO 11
          ASCAT=1001.
          GO TO 80
   11  SS=(D-AD)/(D+AD)
       Q=RR/SS
       SS=AMAX1(0.1,SS)
       Q=AMIN1(AMAX1(0.1,Q),10.)
       Z0=(D-AD)*(D+AD)*TH*0.25/D
       ET=(ETQ*EXP(-AMIN1(1.7,Z0/8.0E3)**6)+1.)*Z0/1.7556E3
       ETT=AMAX1(ET,1.)
       H0=(H0F(R1,ETT)+H0F(R2,ETT))*0.5
       H0=H0+AMIN1(H0,(1.38-ALOG(ETT))*ALOG(SS)*ALOG(Q)*0.49)
       H0=DIM(H0,0.)
       IF(ET .LT. 1.) H0=ET*H0+(1.-ET)*4.343*ALOG(((1.+1.4142/R1)*
     X    (1.+1.4142/R2))**2*(R1+R2)/(R1+R2+2.8284))
       IF(H0 .LE. 15. .OR. H0S .LT. 0.) GO TO 13
  12      H0=H0S
  13   H08=H0
       TH=THA+D*GME
       ASCAT=AHD(TH*D)+4.343*ALOG(47.7*WN*TH**4)-
     X    0.1*(ENS-301.)*EXP(-TH*O/40E3)+H0
  80   RETURN
      END

      FUNCTION AKNFE(V2)
C        KNIFE-EDGE DIFFRACTION
C        THE FRESNEL INTEGRAL AS A FUNCTION OF V**2
C
       IF(V2 .GT. 5.76) GO TO 2
  1       AKNFE=6.02+9.11*SQRT(V2)-1.27*V2
        GO TO 8
  2       AKNFE=12.953+4.343*ALOG(V2)
  8    RETURN
      END

      FUNCTION FHT(X,PK)
C        THE HEIGHT GAIN OVER A SMOOTH SPHERICAL EARTH
C        TO BE USED IN THE *THREE RADII* METH0D
C
       IF(X .LT. 200.) GO TO 2
          FHT=0.05751*X-4.343*ALOG(X)
        IF(X .GE. 2000.) GO TO 8
          W=0.0134*X*EXP(-0.005*X)
          FHT=(1.-W)*FHT+W*(17.372*ALOG(X)-117.)
        GO TO 8
   2      IF(PK .GT. 1.E-5) GO TO 3
          IF(X .GT. 1.) GO TO 4
          FHT=-117.
        GO TO 8
   3      W=-ALOG(PK)
          IF(X*W**3 .GT. 5495.) GO TO 4
          FHT=2.5E-5*X**2/PK-8.686*W-15.
        GO TO 8
   4      FHT=17.372*ALOG(X)-117.
   8   RETURN
      END

      FUNCTION H0F(R,ET)
C       THE H0 FUNCTION FOR SCATTER FIELDS
C
      DIMENSION A(5),B(5)
C
      DATA A( 1 ) ,A( 2) ,A( 3) ,A( 4) ,A( 5)
     X    / 25., 80.,177.,395.,705./
      DATA B(1),B(2),B(3),B(4),B(5)
     X    / 24., 45., 68., 80.,105./
C
       IT=ET
       IF(IT) 2,2,1
  1    IF(IT-5) 5,4,3
  2       IT=1
        GO TO 4
  3       IT=5
  4       Q=0.
        GO TO 6
  5       Q=ET-FLOAT(IT)
  6    X=(1./R)**2
       H0F=4.343*ALOG((A(IT)*X+B(IT))*X+1.)
       IF(Q .NE. 0.)
     X    H0F=(1.-Q)*H0F+Q*4.343*ALOG((A(IT+1)*X+B(IT+1))*X+1.)
       RETURN
      END

      FUNCTION AHD(TD)
C        THE F(TH*D) FUNCTION FOR SCATTER FIELDS
C
      DIMENSION A(3),B(3),C(3)
C
      DATA A(1),A(2),A(3)/133.4,104.6,71.8/
      DATA B(1),B(2),B(3)/0.332E-3,0.212E-3,0.157E-3/
      DATA C(1),C(2),C(3)/-4.343,-1.086,2.171/
C
          I=1
        IF(TD .LE. 10E3) GO TO 1
          I=2
        IF(TD .LE. 70E3) GO TO 1
          I=3
   1   AHD=A(I)+B(I)*TD+C(I)*ALOG(TD)
       RETURN
      END

      FUNCTION AVAR(ZZT,ZZL,ZZC)
C        QUANTILES OF ATTENUATION RELATIVE TO FREE SPACE
C        INCLUDES LONG-TERM TIME VARIABILITY, LOCATION VARIABILITY,
C          AND SITE VARIABILITY
C
      COMMON/PROP/KWX,AREF,MDP,DIST,HG(2),WN,DH,ENS,GME,ZGND,
     X   HE(2),DL(2),THE(2)
         COMPLEX ZGND
      COMMON/PROPV/LVAR,SGC,MDVAR,KLIM
C
      COMMON/SAVE/SAVE(13),KDV,WL,WS,DEXA,DE,VMD,VS0,SGL,
     X   SGTM,SGTP,SGTD,TGTD,GM,GP,CV1,CV2,YV1,YV2,YV3,
     X   CSM1,CSM2,YSM1,YSM2,YSM3,CSP1,CSP2,YSP1,YSP2,YSP3,CSD1,ZD,
     X   CFM1,CFM2,CFM3,CFP1,CFP2,CFP3
C
      DIMENSION BV1(7),BV2(7),XV1(7),XV2(7),XV3(7)
      DIMENSION BSM1(7),BSM2(7),XSM1(7),XSM2(7),XSM3(7)
      DIMENSION BSP1(7),BSP2(7),XSP1(7),XSP2(7),XSP3(7)
      DIMENSION BSD1(7),BZD1(7)
      DIMENSION BFM1(7),BFM2(7),BFM3(7),BFP1(7),BFP2(7),BFP3(7)
C
      LOGICAL WS,WL
C
C          EQUATOR, CON SUB, MAR SUB, DESERT, CON TMP, MAR LND, MAR SEA
C
      DATA  BV1(1), BV1(2), BV1(3), BV1(4), BV1(5), BV1(6), BV1(7)
     X    /  -9.67,  -0.62,   1.26,  -9.21,  -0.62,  -0.39,   3.15/
      DATA  BV2(1), BV2(2), BV2(3), BV2(4), BV2(5), BV2(6), BV2(7)
     X    /   12.7,   9.19,   15.5,   9.05,   9.19,   2.86,  857.9/
      DATA  XV1(1), XV1(2), XV1(3), XV1(4), XV1(5), XV1(6), XV1(7)
     X    /144.9E3,228.9E3,262.6E3, 84.1E3,228.9E3,141.7E3,2222.E3/
      DATA  XV2(1), XV2(2), XV2(3), XV2(4), XV2(5), XV2(6), XV2(7)
     X    /190.3E3,205.2E3,185.2E3,101.1E3,205.2E3,315.9E3,164.8E3/
      DATA  XV3(1), XV3(2), XV3(3), XV3(4), XV3(5), XV3(6), XV3(7)
     X    /133.8E3,143.6E3, 99.8E3, 98.6E3,143.6E3,167.4E3,116.3E3/
      DATA BSM1(1),BSM1(2),BSM1(3),BSM1(4),BSM1(5),BSM1(6),BSM1(7)
     X    / 2.13,     2.66,   6.11,   1.98,   2.68,   6.86,   8.51/
      DATA BSM2(1),BSM2(2),BSM2(3),BSM2(4),BSM2(5),BSM2(6),BSM2(7)
     X    / 159.5,    7.67,   6.65,  13.11,   7.16,  10.38,  169.8/
      DATA XSM1(1),XSM1(2),XSM1(3),XSM1(4),XSM1(5),XSM1(6),XSM1(7)
     X    /762.2E3,100.4E3,138.2E3,139.1E3, 93.7E3,187.8E3,609.8E3/
      DATA XSM2(1),XSM2(2),XSM2(3),XSM2(4),XSM2(5),XSM2(6),XSM2(7)
     X    /123.6E3,172.5E3,242.2E3,132.7E3,186.8E3,169.6E3,119.9E3/
      DATA XSM3(1),XSM3(2),XSM3(3),XSM3(4),XSM3(5),XSM3(6),XSM3(7)
     X    / 94.5E3,136.4E3,178.6E3,193.5E3,133.5E3,108.9E3,106.6E3/
      DATA BSP1(1),BSP1(2),BSP1(3),BSP1(4),BSP1(5),BSP1(6),BSP1(7)
     X    /   2.11,   6.87,  10.08,   3.68,   4.75,   8.58,   8.43/
      DATA BSP2(1),BSP2(2),BSP2(3),BSP2(4),BSP2(5),BSP2(6),BSP2(7)
     X    /  102.3,  15.53,   9.60,  159.3,   8.12,  13.97,   8.19/
      DATA XSP1(1),XSP1(2),XSP1(3),XSP1(4),XSP1(5),XSP1(6),XSP1(7)
     X    /636.9E3,138.7E3,165.3E3,464.4E3, 93.2E3,216.0E3,136.2E3/
      DATA XSP2(1),XSP2(2),XSP2(3),XSP2(4),XSP2(5),XSP2(6),XSP2(7)
     X    /134.8E3,143.7E3,225.7E3, 93.1E3,135.9E3,152.0E3,188.5E3/
      DATA XSP3(1),XSP3(2),XSP3(3),XSP3(4),XSP3(5),XSP3(6),XSP3(7)
     X    / 95.6E3, 98.6E3,129.7E3, 94.2E3,113.4E3,122.7E3,122.9E3/
      DATA BSD1(1),BSD1(2),BSD1(3),BSD1(4),BSD1(5),BSD1(6),BSD1(7)
     X    /  1.224,  0.801,  1.380,  1.000,  1.224,  1.518,  1.518/
      DATA BZD1(1),BZD1(2),BZD1(3),BZD1(4),BZD1(5),BZD1(6),BZD1(7)
     X    /  1.282,  2.161,  1.282,    20.,  1.282,  1.282,  1.282/
      DATA BFM1(1),BFM1(2),BFM1(3),BFM1(4),BFM1(5),BFM1(6),BFM1(7)
     X    /     1.,     1.,     1.,     1.,   0.92,     1.,     1./
      DATA BFM2(1),BFM2(2),BFM2(3),BFM2(4),BFM2(5),BFM2(6),BFM2(7)
     X    /     0.,     0.,     0.,     0.,  0 .25,     0.,     0./
      DATA BFM3(1),BFM3(2),BFM3(3),BFM3(4),BFM3(5),BFM3(6),BFM3(7)
     X    /     0.,     0.,     0.,     0.,   1.77,     0.,     0./
      DATA BFP1(1),BFP1(2),BFP1(3),BFP1(4),BFP1(5),BFP1(6),BFP1(7)
     X    /     1.,   0.93,     1.,   0.93,   0.93,     1.,     1./
      DATA BFP2(1),BFP2(2),BFP2(3),BFP2(4),BFP2(5),BFP2(6),BFP2(7)
     X    /     0.,   0.31,     0.,   0.19,   0.31,     0.,     0./
      DATA BFP3(1),BFP3(2),BFP3(3),BFP3(4),BFP3(5),BFP3(6),BFP3(7)
     X    /     0.,   2.00,     0.,   1.79,   2.00,     0.,     0./
C
      DATA RT,RL/7.8,24./
C
      DATA THIRD/0.3333333/
C
      CURV(C1,C2,X1,X2,X3)=(C1+C2/(1.+((DE-X2)/X3)**2))*
     X   ((DE/X1)**2)/(1.+((DE/X1)**2))
C
      IF(LVAR .EQ. 0) GO TO 60
      IF(LVAR .LT. 5)
     X   GO TO (10,20,30,40),LVAR
C       CLIMATE
      IF(KLIM .GT. 0. AND. KLIM .LE. 7) GO TO 51
         KLIM=5
         KWX=MAX0(KWX,2)
 51   CV1=BV1(KLIM)
      CV2=BV2(KLIM)
      YV1=XV1(KLIM)
      YV2=XV2(KLIM)
      YV3=XV3(KLIM)
      CSM1=BSM1(KLIM)
      CSM2=BSM2(KLIM)
      YSM1=XSM1(KLIM)
      YSM2=XSM2(KLIM)
      YSM3=XSM3(KLIM)
      CSP1=BSP1(KLIM)
      CSP2=BSP2(KLIM)
      YSP1=XSP1(KLIM)
      YSP2=XSP2(KLIM)
      YSP3=XSP3(KLIM)
      CSD1=BSD1(KLIM)
      ZD=BZD1(KLIM)
      CFM1=BFM1(KLIM)
      CFM2=BFM2(KLIM)
      CFM3=BFM3(KLIM)
      CFP1=BFP1(KLIM)
      CFP2=BFP2(KLIM)
      CFP3=BFP3(KLIM)
C       MODE OF VARIABILITY
  40  KDV=MDVAR
      WS=KDV .GE. 20
      IF(WS) KDV=KDV-20
      WL=KDV .GE. 10
      IF(WL) KDV=KDV-10
      IF(KDV .GE. 0 .AND. KDV .LE. 3) GO TO 41
         KDV=0
         KWX=MAX0(KWX,2)
  41  KDV=KDV+ 1
C       FREQUENCY
  30  Q=ALOG(0.133*WN)
      GM=CFM1+CFM2/((CFM3*Q)**2+1.)
      GP=CFP1+CFP2/((CFP3*Q)**2+1.)
C       SYSTEM PARAMETERS
  20  DEXA=SQRT(18E6*HE(1))+SQRT(18E6*HE(2))+(575.7E12/WN)**THIRD
C       DISTANCE
  10  IF(DIST .GE. DEXA) GO TO 11
         DE=130E3*DIST/DEXA
       GO TO 12
  11     DE=130E3+DIST-DEXA
  12  VMD=CURV( CV1 , CV2, YV1, YV2, YV3)
      SGTM=CURV(CSM1,CSM2,YSM1,YSM2,YSM3)*GM
      SGTP=CURV(CSP1,CSP2,YSP1,YSP2,YSP3)*GP
      SGTD=SGTP*CSD1
      TGTD=(SGTP-SGTD)*ZD
      SGL=0.
      VS0=0.
      IF(WL) GO TO 13
         Q=(1.-0.8*EXP(-AMIN1(20.,DIST/50E3)))*DH*WN
         SGL=10.*Q/(Q+13.)
  13  IF(WS) GO TO 14
      VS0=(5.+3.*EXP(-AMIN1(20.,DE/100E3)))**2
  14  CONTINUE
C
      LVAR=0
C
  60  CONTINUE
      ZT=ZZT
      ZL=ZZL
      ZC=ZZC
      GO TO (600,601,602,603),KDV
  600    ZT=ZC
  601    ZL=ZC
       GO TO 603
  602 ZL=ZT
  603 CONTINUE
      IF(ABS(ZT) .GT. 3.10) GO TO 605
      IF(ABS(ZL) .GT. 3.10) GO TO 605
      IF(ABS(ZC) .GT. 3.10) GO TO 605
      GO TO 608
  605    KWX=MAX0(KWX,1)
  608 CONTINUE
      IF(ZT .GT. 0.) GO TO 611
         SGT=SGTM
       GO TO 618
  611    IF(ZT .GT. ZD) GO TO 612
         SGT=SGTP
       GO TO 618
  612    SGT=SGTD+TGTD/ZT
  618 CONTINUE
      VS=VS0+(SGT*ZT)**2/(RT+ZC**2)+(SGL*ZL)**2/(RL+ZC**2)
      GO TO (620,621,622,623),KDV
  620    YR=0.
         SGC=SQRT(SGT**2+SGL**2+VS)
       GO TO 628
  621    YR=SGT*ZT
         SGC=SQRT(SGL**2+VS)
       GO TO 628
  622    YR=SQRT(SGT**2+SGL**2)*ZT
         SGC=SQRT(VS)
       GO TO 628
  623    YR=SGT*ZT+SGL*ZL
         SGC=SQRT(VS)
  628 CONTINUE
C
      AVAR=AREF-VMD-YR-SGC*ZC
      IF(AVAR .LT. 0.) AVAR=AVAR*(29.-AVAR)/(29.-10.*AVAR)
      RETURN
      END

      FUNCTION QERF(Z)
C        THE STANDARD NORMAL COMPLEMENTARY PROBABILITY
C        APPROXIMATION DUE TO C. HASTINGS, JR.
C        MAX ERROR 7.5E-8
C
      DATA B1,B2,B3,B4,B5/0.319381530,-0.356563782,1.781477937,
     X   -1.821255987,1.330274429/
      DATA RP,RRT2PI/4.317008,0.398942280/
C
      X=Z
      T=ABS(X)
      IF(T .LT. 10.) GO TO 1
      QERF=0.
      GO TO 2
  1   T=RP/(T+RP)
      QERF=EXP(-0.5*X**2)*RRT2PI*((((B5*T+B4)*T+B3)*T+B2)*T+B1)*T
  2   IF(X .LT. 0.) QERF=1.-QERF
      RETURN
      END

      FUNCTION QERFI(Q)
C        THE INVERSE OF QERF, GIVES THE STANDARD NORMAL DEVIATE AS A
C           FUNCTION OF THE COMPLEMENTARY PROBABILITY
C        TRUNCATED AT 0.000001 AND 0.999999
C        APPROXIMATION DUE TO C. HASTINGS, JR.
C        MAX ERROR 4.5E-4
C
      DATA C0,C1,C2/2.515516698,0.802853,0.010328/
      DATA D1,D2,D3/1.432788,0.189269,0.001308/
C
      X=0.5-Q
      T=AMAX1(0.5-ABS(X),0.000001)
      T=SQRT(-2.*ALOG(T))
      QERFI=T-((C2*T+C1)*T+C0)/(((D3*T+D2)*T+D1)*T+1.)
      IF(X .LT. 0.) QERFI=-QERFI
      RETURN
      END

      SUBROUTINE QLRPS(FMHZ,ZSYS,EN0,IPOL,EPS,SGM)
C        PREPARES PARAMETERS
C        SETS--
C          WN,ENS,GME,ZGND
C
      COMMON/PROP/KWX,AREF,MDP,DIST,HG(2),WN,DH,ENS,GME,ZGND,
     X   HE(2),DL(2),THE(2)
      COMPLEX ZGND
C
      COMPLEX ZQ
C
      DATA GMA/157E-9/
C
       WN=FMHZ/47.7
       ENS=EN0
       IF(ZSYS .NE. 0.) ENS=ENS*EXP(-ZSYS/9460.)
       GME=GMA*(1.-0.04665*EXP(ENS/179.3))
       ZQ=CMPLX(EPS,376.62*SGM/WN)
       ZGND=CSQRT(ZQ-1.)
       IF(IPOL .NE. 0) ZGND=ZGND/ZQ
       RETURN
      END

      SUBROUTINE QLRA( KST, KLIMX,MDVARX)
         DIMENSION KST(2)
C        PREPARES THE LONGLEY-RICE MODEL IN THE AREA PREDICTION MODE
C
      COMMON/PROP/KWX,AREF,MDP,DIST,HG(2),WN,DH,ENS,GME,ZGND,
     X   HE(2),DL(2),THE(2)
         COMPLEX ZGND
      COMMON/PROPV/LVAR,SGC,MDVAR,KLIM
C
       DO 10 J=1,2
         IF(KST(J)-1) 11,12,13
  11        HE(J)=HG(J)
          GO TO 15
  12        Q=4.
            GO TO 14
  13        Q=9.
  14        IF(HG(J) .LT. 5.) Q=Q*SIN(0.3141593*HG(J))
            HE(J)=HG(J)+(1.+Q)*EXP(-AMIN1(20.,2.*HG(J)/AMAX1(1E-3,DH)))
  15     Q=SQRT(2.*HE(J)/GME)
         DL(J)=Q*EXP(-0.07*SQRT(DH/AMAX1(HE(J),5.)))
         THE(J)=(0.65*DH*(Q/DL(J)-1.)-2.*HE(J))/Q
  10   CONTINUE
C
       MDP=1
       LVAR=MAX0(LVAR,3)
       IF(MDVARX .LT. 0) GO TO 21
          MDVAR=MDVARX
          LVAR=MAX0(LVAR,4)
  21   IF(KLIMX .LE. 0) GO TO 22
          KLIM=KLIMX
          LVAR=5
  22   CONTINUE
       RETURN
      END

      SUBROUTINE QLRPFL(PFL, KLIMX, MDVARX )
         DIMENSION PFL(5)
C
C        SETS UP AND RUNS THE LONGLEY-RICE MODEL IN THE POINT-TO-POINT
C          MODE USING THE TERRAIN PROFILE IN PFL.
C            PFL(1)=ENP, PFL(2)=XI, PFL(3)=Z(0),...
C
      COMMON/PROP/KWX,AREF,MDP,DIST,HG(2),WN,DH,ENS,GME,ZGND,
     X   HE(2),DL(2),THE(2)
         COMPLEX ZGND
      COMMON/PROPV/LVAR,SGC,MDVAR,KLIM
C
      DIMENSION XL(2)
C
       DIST=PFL(1)*PFL(2)
       NP=PFL(1)
       CALL HZNS(PFL)
C        FIND DELTA H
       DO 11 J=1,2
  11   XL(J)=AMIN1(15.*HG(J),0.1*DL(J))
       XL(2)=DIST-XL(2)
       DH=DLTHX(PFL,XL(1),XL(2))
C        FIND EFFECTIVE HEIGHTS HE
       IF(DL(1)+DL(2) .LT. 1.5*DIST) GO TO 25
C        LINE-OF-SIGHT
          CALL ZLSQ1(PFL,XL(1),XL(2),ZA,ZB)
          HE(1)=HG(1)+DIM(PFL(3),ZA)
          HE(2)=HG(2)+DIM(PFL(NP+3),ZB)
          DO 21 J=1,2
  21      DL(J)=SQRT(2.*HE(J)/GME)*EXP(-0.07*SQRT(DH/AMAX1(HE(J),5.)))
          Q=DL(1)+DL(2)
          IF(Q .GT. DIST) GO TO 23
          Q=(DIST/Q)**2
          DO 22 J=1,2
          HE(J)=HE(J)*Q
  22      DL(J)=SQRT(2.*HE(J)/GME)*EXP(-0.07*SQRT(DH/AMAX1(HE(J),5.)))
  23    GO TO 28
C        TRANSHORIZ0N
  25      CALL ZLSQ1(PFL,XL(1),0.9*DL(1),ZA,Q)
          CALL ZLSQ1(PFL,DIST-0.9*DL(2),XL(2),Q,ZB)
          HE(1)=HG(1)+DIM(PFL(3),ZA)
          HE(2)=HG(2)+DIM(PFL(NP+3),ZB)
  28   CONTINUE
C
       MDP=-1
       LVAR=MAX0(LVAR,3)
       IF(MDVARX .LT. 0) GO TO 31
          MDVAR=MDVARX
          LVAR=MAX0(LVAR,4)
  31   IF(KLIMX .LE. 0) GO TO 32
          KLIM=KLIMX
          LVAR=5
  32   CONTINUE
C
       CALL LRPROP(0.)
C
       RETURN
      END

      FUNCTION DLTHX(PFL,X1,X2)
         DIMENSION PFL(5)
C
C        COMPUTES THE TERRAIN IRREGULARITY PARAMETER DH FROM THE
C        PROFILE PFL BETWEEN POINTS AT X1 .LT. X2.
C
      DIMENSION S(247)
C
       NP=PFL(1)
       XA=X1/PFL(2)
       XB=X2/PFL(2)
       DLTHX=0.
       IF(XB-XA .LT. 2.) GO TO 80
       KA=0.1*(XB-XA+8.)
       KA=MIN0(MAX0(4,KA),25)
       N=10*KA-5
       KB=N-KA+1
       SN=N-1
       S(1)=SN
       S(2)=1.
       XB=(XB-XA)/SN
       K=XA+1.
       XA=XA-FLOAT(K)
       DO 10 J=1,N
  11   IF(XA .LE. 0.) GO TO 12
       IF(K .GE. NP) GO TO 12
         XA=XA-1.
         K=K+1
       GO TO 11
  12   S(J+2)=PFL(K+3)+(PFL(K+3)-PFL(K+2))*XA
  10   XA=XA+XB
       CALL ZLSQ1(S,0.,SN,XA,XB)
       XB=(XB-XA)/SN
       DO 15 J=1,N
       S(J+2)=S(J+2)-XA
  15   XA=XA+XB
C
       DLTHX=QTILE(N,S(3),KA)-QTILE(N,S(3),KB)
       DLTHX=DLTHX/(1.-0.8*EXP(-AMIN1(20.,(X2-X1)/50E3)))
  80   RETURN
      END

      SUBROUTINE HZNS(PFL)
         DIMENSION PFL(5)
C
C        TO FIND HORIZ0NS FROM ANTENNAS WITH HEIGHTS HG AT THE TWO
C          ENDS OF THE PROFILE PFL.
C            PFL(1)=ENP, PFL(2)=XI, PFL(3)=Z(D), ...
C        OUTPUT--DISTANCES DL, TAKE-OFF ANGLES THE.
C          DL=DIST IF THE PATH IS LINE OF SIGHT
C
      COMMON/PROP/KWX,AREF,MDP,DIST,HG(2),WN,DH,ENS,GME,ZGND,
     X   HE(2),DL(2),THE(2)
      COMPLEX ZGND
C
      LOGICAL WQ
C
       NP=PFL(1)
       XI=PFL(2)
       ZA=PFL(3)+HG(1)
       ZB=PFL(NP+3)+HG(2)
       QC=0.5*GME
       Q=QC*DIST
       THE(2)=(ZB-ZA)/DIST
       THE(1)=THE(2)-Q
       THE(2)=-THE(2)-Q
       DL(1)=DIST
       DL(2)=DIST
       IF(NP .LT. 2) GO TO 18
       SA=0.
       SB=DIST
       WQ=.TRUE.
       DO 10 I=2,NP
         SA=SA+XI
         SB=SB-XI
         Q=PFL(I+2)-(QC*SA+THE(1))*SA-ZA
         IF(Q .LE. 0.) GO TO 11
            THE(1)=THE(1)+Q/SA
            DL(1)=SA
            WQ=.FALSE.
  11        IF(WQ) GO TO 10
            Q=PFL(I+2)-(QC*SB+THE(2))*SB-ZB
            IF(Q .LE. 0.) GO TO 10
            THE(2)=THE(2)+Q/SB
            DL(2)=SB
  10   CONTINUE
C
  18   RETURN
      END

      FUNCTION QTILE(NN,A,IR)
         DIMENSION A(NN)
C
C        REORDERS A SO THAT A(J),J=1 ... IR ARE ALL .GE.
C          ALL A(I),I=IR ... NN. IN PARTICULAR, A(IR) WILL HAVE THE SAME
C          VALUE IT WOULD HAVE IF A WERE COMPLETELY SORTED IN
C          DESCENDING ORDER.
C        RETURNS QTILE=A(IR)
C
       M=1
       N=NN
       K=MIN0(MAX0(1,IR),N)
  10   CONTINUE
         Q=A(K)
         I0=M
         J1=N
  11     CONTINUE
           DO 12 I=I0,N
           IF(A(I) .LT. Q) GO TO 13
  12       CONTINUE
           I=N
  13       J=J1
           DO 14 JJ=M,J1
           IF(A(J) .GT. Q) GO TO 15
  14       J=J-1
           J=M
  15       IF(I .GE. J) GO TO 16
           R=A(I)
           A(I)=A(J)
           A(J)=R
           I0=I+1
           J1=J-1
           GO TO 11
  16     IF(I .GE. K) GO TO 17
         A(K)=A(I)
         A(I)=Q
         M=I+1
         GO TO 10
  17     IF(J .LE. K) GO TO 20
         A(K)=A(J)
         A(J)=Q
         N=J-1
         GO TO 10
  20   QTILE=Q
       RETURN
      END

      SUBROUTINE ZLSQ1(Z,X1,X2,Z0,ZN)
         DIMENSION Z(5)
C
C        LINEAR LEAST SQUARES FIT BE'IWEEN X1, X2 TO THE FUNCTION
C           DESCRIBED BY Z--
C          Z(1)=EN, NUMBER OF INTERVALS, Z(2) =XI, INTERVAL LENGTH,
C          Z(J+3), J=0,...,EN, FUNCTION VALUES.
C        OUTPUT-- VALUES OF THE LINE, Z0 AT 0, ZN AT XT.
C
       XN=Z(1)
       XA=AINT(DIM(X1/Z(2),0.))
       XB=XN-AINT(DIM(XN,X2/Z(2)))
       IF(XB .GT. XA) GO TO 1
       XA=DIM(XA,1.)
       XB=XN-DIM(XN,XB+1.)
  1    JA=XA
       JB=XB
       N=JB-JA
       XA=XB-XA
       X=-0.5*XA
       XB=XB+X
       A=0.5*(Z(JA+3)+Z(JB+3))
       B=0.5*(Z(JA+3)-Z(JB+3))*X
       IF(N .LT. 2) GO TO 11
       DO 10 I=2,N
         JA=JA+1
         X=X+1.
         A=A+Z(JA+3)
         B=B+Z(JA+3)*X
  10   CONTINUE
  11   A=A/XA
       B=B*12./((XA*XA+2.)*XA)
       Z0=A-B*XB
       ZN=A+B*(XN-XB)
       RETURN
      END

C     PROGRAM QKAREA
C        *QUICK AREA*
C        TO ILLUSTRATE THE USE OF THE LONGLEY-RICE MODEL
C          IN THE AREA PREDICTION MODE
C
C        INPUT IS IN 10-COL FIELDS, THE FIRST OF WHICH IS
C          A SEQUENCE OF DIGITS
C        IN PARTICULAR,
C          COL 1 IS THE *EXECUTE* COLUMN--A NON-ZERO DIGIT
C             WILL FORCE OUTPUT
C          COL 2 INDICATES THE CARD TYPE--
C
C                          COL
C                          12      11,...
C            STOP-         X0         (OR A BLANK CARD)
C            TITLE-        X1         (NEXT CARD HAS 60-COL TITLE)
C            DISTANCES-    X2      D0,D1,DS1,D2,DS2
C            RELIABILITY-  X3V     QT,QL
C            CONFIDENCE-   X4      QC1,QC2,...
C            ENVIRONMENT-  X5C     DH,N0,ZS,EPS,SGM
C            SYSTEM-       X6NPSS  FMHZ,HG1,HG2
C            (ALTERNATE)   X7NPSS  FMHZ,HG1,HG2,DH,NS,EPS,SGM
C            EXECUTE-      X8
C            RESET-        X9
C

      COMMON/PROP/KWX,AREF,MDP,DIST,HG(2),WN,DH,ENS,GME,ZGND,
     X   HE(2),DL(2),THE(2)
         COMPLEX ZGND
      COMMON/PROPV/LVAR,SGC,MDVAR,KLIM
C
      COMMON/PROPA/DLSA,DX,AEL,AK1,AK2,AED,EMD,AES,EMS,DLS(2),DLA,THA
      COMMON/SAVE/SAVE(50)
C
      DIMENSION JIN(6),XIN(7)
      DIMENSION ITL(15)
      DIMENSION KST(2)
      DIMENSION QC(7),ZC(7),XLB(7)
C
      LOGICAL WQIT,WCON,WTL
C
C        THE I/O UNITS ARE DEFINED HERE
      DATA KIN,KOT/5,6/
C
      DATA GMA/157E-9/
C
      DATA DB/8.685890/
      DATA AKM/1000./
C
      WQIT=.FALSE.
      WCON=.TRUE.
      GO TO 190
C
  10  CONTINUE
C       READ INPUT SEQUENCE
C
1000  FORMAT(6I1,4X,7F10.0)
1001  FORMAT(15A4)
C
      JIN(1)=0
      JIN(2)=0
      READ(KIN,1000) JIN,XIN
      WCON=JIN(1) .EQ. 0
      JQ=JIN(2)
      IF(JQ .NE. 0)
     X   GO TO (110,120,130,140,150,160,170,180,190),JQ
C
         WQIT=.TRUE.
      GO TO 20
 110    CONTINUE
        READ(KIN,1001) ITL
        WTL=.TRUE.
      GO TO 20
 120    CONTINUE
        XIN(1)=DIM(XIN(1),0.)
        Q=XIN(2)-XIN(1)
        IF(Q .GT. 0.) GO TO 121
           IF(XIN(1) .EQ. 0.) GO TO 128
           D0=XIN(1)
           DS=0.
           DSC=0.
           ND=1
           NDC=0
      GO TO 128
 121    IF(XIN(3) .LE. 0.) XIN(3)=AMAX1(1.,AINT(Q/20.+0.5))
        IF(XIN(1) .LE. 0.) XIN(1)=XIN(3)
        D0=XIN(1)
        DS=XIN(3)
        DSC=DS
        ND=DIM(XIN(2),XIN(1))/DS+1.75
        NDC=0
      IF(XIN(4) .LE. XIN(2)) GO TO 128
        IF(XIN(5) .LE. 0.) XIN(5)=5.*XIN(3)
        DSC=XIN(5)
        JQ=(XIN(4)-XIN(2))/DSC+0.75
        NDC=ND
        ND=ND+JQ
 128  GO TO 20
 130    CONTINUE
        MDVAR=MIN0(JIN(3),3)
        LVAR=MAX0(LVAR,4)
        QT=50.
        QL=50.
        ZT=0.
        ZL=0.
        IF(XIN(1) .LE. 0.) GO TO 131
           QT=XIN(1)
           ZT=QERFI(QT/100.)
 131       IF(XIN(2) .LE. 0.) GO TO 138
           QL=XIN(2)
           ZL=QERFI(QL/100.)
 138  GO TO 20
 140    CONTINUE
        NC=0
        DO 141 JC=1,7
          IF(XIN(JC) .LE. 0.) GO TO 141
          NC=NC+1
          QC(NC)=XIN(JC)
          ZC(NC)=QERFI(QC(NC)/100.)
 141    CONTINUE
        IF(NC .GT. 0) GO TO 148
           NC=1
           QC(1)=50.
           ZC(1)=0.
 148  GO TO 20
 150    CONTINUE
        IF(JIN(3) .LE. 0) GO TO 151
           KLIM=JIN(3)
           LVAR=5
 151    IF(XIN(1) .GE. 0.) DH=XIN(1)
        IF(XIN(2) .LE. 0.) GO TO 152
           EN0=XIN(2)
           ZSYS=XIN(3)
 152    IF(XIN(4) .LE. 0.) GO TO 158
           EPS=XIN(4)
           SGM=XIN(5)
 158  GO TO 20
 160    CONTINUE
        IF(JIN(3) .EQ. 1) GO TO 161
           IPOL=MIN0(JIN(4),1)
           KST(1)=MIN0(JIN(5),2)
           KST(2)=MIN0(JIN(6),2)
 161    IF(XIN( 1) .GT. 0.) FMHZ=XIN( 1)
        IF(XIN(2) .GT. 0.) HG(1)=XIN(2)
        IF(XIN(3) .GT. 0.) HG(2)=XIN(3)
      GO TO 20
 170  CONTINUE
      IF(JIN(3) .EQ. 1) GO TO 171
         IPOL=MIN0(JIN(4),1)
         KST(1)=MIN0(JIN(5),2)
         KST(2)=MIN0(JIN(6),2)
 171  IF(XIN(1) .GT. 0.) FMHZ=XIN(1)
      IF(XIN(2) .GT. 0.) HG(1)=XIN(2)
      IF(XIN(3) .GT. 0.) HG(2)=XIN(3)
      IF(XIN(4) .GE. 0.) DH=XIN(4)
      IF(XIN(5) .LE. 0.) GO TO 172
        EN0=XIN(5)
        ZSYS=0.
 172  IF(XIN(6) .LE. 0.) GO TO 178
        EPS=XIN(6)
        SGM=XIN(7)
 178  GO TO 20
 180    CONTINUE
        WCON=.FALSE.
      GO TO 20
 190    CONTINUE
        FMHZ=100.
        HG(1)=3.
        HG(2)=3.
        DH=90.
        EN0=301.
        ZSYS=0.
        EPS=15.
        SGM=0.005
        IPOL=1
        KST(1)=0
        KST(2)=0
        KLIM=5
        MDVAR=3
        LVAR=5
        NC=3
        QC(1)=50.
        QC(2)=90.
        QC(3)=10.
        QT=50.
        QL=50.
        ZC(1)=0.
        ZC(2)=-1.28155
        ZC(3)= 1.28155
        ZT=0.
        ZL=0.
        D0=10.
        DS=10.
        DSC=50.
        ND=22
        NDC=15
        WTL=.FALSE.
C
  20  CONTINUE
      IF(WCON) GO TO 30
C
C       EXECUTION
C
      KWX=0
      CALL QLRPS(FMHZ,ZSYS,EN0,IPOL,EPS,SGM)
      CALL QLRA(KST,-1,-1)
C
C     WRITE HEADING
2001  FORMAT(1H1/1H0)
2002  FORMAT(1H )
2010  FORMAT(3X,
     .62HAREA PREDICTIONS FROM THE LONGLEY-RICE MODEL, VERSION 1.2.1   )
2011  FORMAT(3X,15A4)
2015  FORMAT(12X,9HFREQUENCY,F12.0,4H MHZ)
2016  FORMAT(6X,15HANTENNA HEIGHTS,2F8.1,2H M)
2017  FORMAT(4X,17HEFFECTIVE HEIGHTS,2F8.1,
     .   12H M  (SITING=,I1,1H,,I1,1H))
2018  FORMAT(5X,16HTERRAIN, DELTA H,F12.0,2H M)
C
c      WRITE(KOT,2001)
       IF(WTL) GO TO 211
       WRITE(KOT,2010)
       GO TO 212
  211  WRITE(KOT,2011) ITL
  212  WRITE(KOT,2002)
       WRITE(KOT,2002)
       WRITE(KOT,2015) FMHZ
       WRITE(KOT,2016) HG
       WRITE(KOT,2017) HE,KST
       WRITE(KOT,2018) DH
       WRITE(KOT,2002)
C
2021   FORMAT(3X,4HPOL=,I1,6H, EPS=,F3.0,6H, SGM=,F6.3,4H S/M)
2022   FORMAT(3X,5HCLIM=,I1,5H, N0=,F4.0,5H, NS=,F4.0,4H, K=,F6.3)
C
       Q=GMA/GME
       WRITE(KOT,2021) IPOL,EPS,SGM
       WRITE(KOT,2022) KLIM,EN0,ENS,Q
       WRITE(KOT,2002)
C
2030  FORMAT(3X,22HSINGLE-MESSAGE SERVICE)
2031  FORMAT(3X,18HACCIDENTAL SERVICE/
     .   8X,F5.1,27H PER CENT TIME AVAILABILITY)
2032  FORMAT(3X,14HMOBILE SERVICE/
     .   8X,21HREQUIRED RELIABILITY-,F5.1,9H PER CENT)
2033  FORMAT(3X,17HBROADCAST SERVICE/
     .   8X,21HREQUIRED RELIABILITY-,F5.1,14H PER CENT TIME/
     .   29X,F5.1,19H PER CENT LOCATIONS)
C
       IF(MDVAR .NE. 0)
     X   GO TO (231,232,233),MDVAR
C
         WRITE(KOT,2030)
       GO TO 238
  231    WRITE(KOT,2031) QT
       GO TO 238
  232    WRITE(KOT,2032) QT
       GO TO 238
  233    WRITE(KOT,2033) QT,QL
  238 WRITE(KOT,2002)
C
2040  FORMAT(3X,
     .62HESTIMATED QUANTILES OF BASIC TRANSMISSION LOSS(DB)            )
2041  FORMAT(7X, 4HDIST, 5X , 4HFREE,4X, 15HWITH CONFIDENCE/
     .   8X,2HKM,5X,6HSPACE ,7F8.1)
2045  FORMAT(2X,3F9.1,6F8.1)
C
C       COMPUTE AND PRINT VALUES
      WRITE(KOT,2040)
      WRITE(KOT,2002)
      WRITE(KOT,2041) (QC(JC),JC=1,NC)
      WRITE(KOT,2002)
      DT=DS
      D=D0
      DO 240 JD=1,ND
        LVAR=MAX0(1,LVAR)
        CALL LRPROP(D*AKM)
        FS=DB*ALOG(2.*WN*DIST)
        DO 241 JC=1,NC
 241    XLB(JC)=FS+AVAR(ZT,ZL,ZC(JC))
        WRITE(KOT,2045) D,FS,(XLB(JC),JC=1,NC)
        IF(JD .EQ. NDC) DT=DSC
        D=D+DT
 240  CONTINUE
C
2081  FORMAT(3X,
     .62H**WARNING- SOME PARAMETERS ARE NEARLY OUT OF RANGE.           /
     .   3X,
     .62H RESULTS SHOULD BE USED WITH CAUTION.                         )
2082  FORMAT(3X,
     .62H**NOTE- DEFAULT PARAMETERS HAVE BEEN SUBSTITUTED              /
     .   3X,
     .62H FOR IMPOSSIBLE ONES.                                         )
2083  FORMAT(3X,
     .62H**WARNING- A COMBINATION OF PARAMETERS IS OUT OF RANGE.       /
     .   3X,
     .62H RESULTS ARE PROBABLY INVALID.                                )
2084  FORMAT(3X,
     .62H**WARNING- SOME PARAMETERS ARE OUT OF RANGE.                  /
     .   3X,
     .62H RESULTS ARE PROBABLY INVALID.                                )
C
       IF(KWX .EQ. 0) GO TO 28
C
C        PRINT ERROR MESSAGES
          WRITE(KOT,2002)
          GO TO (281,282,283,284),KWX
 281      WRITE(KOT,2081)
        GO TO 28
 282      WRITE(KOT,2082)
        GO TO 28
 283      WRITE(KOT,2083)
        GO TO 28
 284      WRITE(KOT,2084)
 28    CONTINUE
C
 30    CONTINUE
       IF(.NOT. WQIT) GO TO 10
C
       STOP
      END
