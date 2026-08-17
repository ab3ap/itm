#!/usr/bin/env python

# This is the Fortran program qkpfl, QuicK ProFiLe, distributed with the
# Longley Rice ITM v1.2.2, which I converted to python changing as little as
# possible.  The hardest part was implementing the equivalent of fortran
# FORMAT statements.  I underestimated how powerful FORMATs are!
#
# To run it, an input file is needed as described in the original comments
# below.  Use the file demoPflData included with this distribution:
#
#   demoPfl < demoPflData
#
# to get a sample run.  The code shows how you would use it with your code.
#
# Mike Markowski, mike.ab3ap@gmail.com
# Dec 2022

from scipy.stats import norm
import itm
import numpy as np
import sys

#     PROGRAM QKPFL
#        *QUICK PROFILE*
#        TO ILLUSTRATE THE USE OF THE LONGLEY-RICE MODEL
#          IN THE POINT-TO-POINT MODE
#
#        INPUT IS IN 10-COL FIELDS, THE FIRST OF WHICH IS
#          A SEQUENCE OF DIGITS
#        IN PARTICULAR,
#          COL 1 IS THE *EXECUTE* COLUMN--A NON-ZERO DIGIT
#             WILL FORCE OUTPUT
#          COL 2 INDICATES THE CARD TYPE--
#                          COL
#                          12      11,..
#            STOP-         X0         (OR A BLANK CARD)
#            TITLE-        X1         (NEXT CARD HAS 60-COL TITLE)
#            PROFILE-      X2T     D,XI,ZSC
#               -WITH PROFILE ELEVATIONS ON THE FOLLOWING CARDS
#            RELIABILITY-  X3      QR1,QR2,..
#            CONFIDENCE-   X4      QC1,QC2,..
#            PARAMETERS-   X7PC    FMHZ,HG1,HG2,N0,NS,EPS,SGM
#            EXECUTE-      X8
#            RESET-        X9
#
#        PROFILE ELEVATION CARDS HAVE 5-COL FIELDS, THE FIRST OF WHICH
#           IS A SEQUENCE OF DIGITS, ALL IN THE FORMAT
#                          ENN  P0,P1,..
#

qc = np.zeros(7)
qr = np.zeros(7)
zc = np.zeros(7)
zr = np.zeros(7)

# Create empty dictionaries.
itm.dictionaries()
itm.mdvarSet(itm.mRandom)
itm.mdvarSet(itm.mPfl)

mzpfl = 600
gammaA = 1/itm.rEarth_m
wqit = False
wcon = True
f_MHz = 100
hg = np.array([3., 3.])
en0 = 310
ens0 = 0
eps = 15
sgm = 0.005
ipol = 1
klim = itm.cContTemp # Continental temperate.
nc = nr = 3
qc[:3] = np.array([50, 90, 10])
qr[:3] = np.array([50, 90, 10])
zc[:3] = norm.isf(qc[:3]/100) # 0., -1.28155, 1.28155
zr[:3] = norm.isf(qr[:3]/100) # 0., -1.28155, 1.28155
dkm = xkm = 0
npt = -1
wpf = False
wtl = False
pfl = np.zeros(mzpfl+3)

while True:
#        READ INPUT SEQUENCE
#
#1000   FORMAT(6I1,4X,7F10.0)
#       READ(KIN,1000) JIN,XIN
    # Replicate Fortran 6I1 format.
    jin = np.zeros(6, dtype=int)
    try:
        line = input()
    except EOFError:
        sys.exit()
    for i in range(min(6, len(line))):
        jin[i] = int(line[i]) if line[i].isnumeric() else 0
    line = line[10:] # Chop off 6I1,4X.
    wcon = (jin[0] == 0) # Is 1st column zero.
    jq = jin[1] # 2nd column of input line.

    # Replicate Fortran 7F10.0 format.
    n = int(np.ceil(len(line)/10)) # Ten or fewer F10.0 floats follow.
    if jq > 0:
        xin = np.zeros(7)
        for i in range(n):
            xin[i] = float(line[:10].replace(' ', '0'))
            line = line[10:]

    if jq == 0: # Quit program after printing info.
        wqit = True
    elif jq == 1: # Read test title.
        itl = input()
        wtl = True # Flag, read in test title.
    elif jq == 2:
        wptl = (jin[2] != 0)
        if wptl: # Read terrain profile title.
            iptl = input()
        zsc = xin[2]
        if zsc <= 0:
            zsc = 1
        npt = -1
        wpf = True
        while wpf:
#1002       Emulate: FORMAT(I1,I2,2X,15F5.0)
            line = input()
            jq = 0 if line[0]==' ' else int(line[0])      # I1
            nq = 0 if line[1:3]=='  ' else int(line[1:3]) # I2
            line = line[5:]                               # Chop I1,I2,2X
            n = int(np.ceil(len(line)/5))                 # Float count.
            qfl = np.zeros(min(15, n))
            for i in range(n):                            # 15F5.0
                qfl[i] = float(line[:5])
                line = line[5:]

            wpf = (jq == 0 )
            nq = min(nq,15)
            if nq <= 0:
                wpf = False
            else:
                for jq in range(nq):
                    npt += 1
                    if npt <= mzpfl:
                        pfl[npt+2] = qfl[jq]*zsc
        pfl[0] = npt
        dkm = xin[0]
        xkm = xin[1]
        wpf = (npt > 0)
        if wpf:
            if dkm <= 0:
                dkm = xkm*pfl[0]
            if xkm <= 0:
                xkm = dkm/pfl[0]
            pfl[1] = dkm*1e3/pfl[0]
            wpf = ((npt <= mzpfl) and (dkm > 0)
                and (abs(dkm-xkm*pfl[0]) < 0.5*xkm))
        pfl = pfl[:3+npt]
    elif jq == 3:
        nr = 0
        for jr in range(7):
            if xin[jr] > 0:
                qr[nr] = xin[jr]
                zr[nr] = norm.isf(qr[nr]/100)
                nr += 1
        if nr <= 0:
            nr = 0
            qr[0] = 50
            zr[0] = 0
    elif jq == 4:
        nc = 0
        for jc in range(7):
            if xin[jc] > 0:
                qc[nc] = xin[jc]
                zc[nc] = norm.isf(qc[nc]/100)
                nc += 1
        if nc <= 0:
            nc = 0
            qc[0] = 50
            zc[0] = 0
    elif jq == 5:
        pass
    elif jq == 6:
        pass
    elif jq == 7:
        ipol = min(jin[2], 1)
        if jin[3] > 0:
            klim = min(jin[3], 7)
            match jin[3]:
                case 1: klim = itm.cEquatorial
                case 2: klim = itm.cContSubtrop
                case 3: klim = itm.cMariSubtrop
                case 4: klim = itm.cDesert
                case 5: klim = itm.cContTemp
                case 6: klim = itm.cMariTempLand
                case 7: klim = itm.cMariTempSea
                case _: klim = itm.cContTemp
        if xin[0] > 0:
            f_MHz = xin[0]
        if xin[1] > 0:
            hg[0] = xin[1]
        if xin[2] > 0:
            hg[1] = xin[2]
        if xin[4] > 0:
            ens0 = xin[4]
        if xin[3] > 0:
            en0 = xin[3]
            ens0 = 0
        if xin[5] > 0:
            eps = xin(6)
            sgm = xin(7)
    elif jq == 8:
        wcon = False
    elif jq == 9:
        f_MHz = 100
        hg[0] = hg[1] = 3.
        en0 = 310
        ens0 = 0
        eps = 15
        sgm = 0.005
        ipol = 1
        klim = itm.cContTemp # Continental temperate.
        nc = nr = 3
        qc[:3] = np.array([50, 90, 10])
        qr[:3] = np.array([50, 90, 10])
        zc[:3] = norm.isf(qc[:3]/100) # 0., -1.28155, 1.28155
        zr[:3] = norm.isf(qr[:3]/100) # 0., -1.28155, 1.28155
        dkm = xkm = 0
        npt = -1
        wpf = False
        wtl = False

    if wcon:
        continue
#
#           EXECUTE
#
#
#           WRITE HEADING

    print('\n')
    if wtl: # Read in test title.
        print('   %s' % itl)
    else:
        print('LINK PREDICTIONS FROM THE LONGLEY-RICE MODEL, VERSION 1.2.2')
    print('')
    if wptl: # Read in profile title.
        print('   %s' % iptl)
        print('')
    print('')
    print('%21s%12.1f KM' % ('DISTANCE', dkm))
    print('%21s%12.1f MHZ' % ('FREQUENCY', f_MHz))
    print('%21s%8.1f%8.1f M' % ('ANTENNA HEIGHTS', hg[0], hg[1]))
    if not wpf:
        print('')
        print('   PROFILE- NP=%4d, XI=%f6.3 KM' % (npt, xkm))
        print('      ERRORS--ANALYSIS CANNOT CONTINUE')
        continue

    itm.prop['hg'] = hg        # Height of antenna structure in m.
    itm.propv['klim'] = klim
    #   C a l c u l a t e   C o n s t a n t s
    zsys = 0
    q = ens0
    if q <= 0:
#       ja = 3 + 0.1*pfl[0]
#       jb = npt - ja + 6 # XXX '+6' seems wrong.
        ja = 2 + 0.1*pfl[0]
        jb = 2 + 0.9*pfl[0]
        zsys = pfl[ja:jb+1].mean()
        q = en0
    pol = 'h' if ipol==0 else 'v'
    itm.lrPrep(f_MHz, zsys, q, pol, eps, sgm)

    #   P e r f o r m   P r o p a g a t i o n   C a l c s
    n = int(pfl[0]) + 1
    rng = pfl[1]*np.arange(n) # Arbitrary distances are also allowed.
    itm.lrProfile(rng, pfl[2:2+n])
#   el1 = np.degrees(itm.prop['the'][0]) # Needed with antenna patterns.
#   el2 = np.degrees(itm.prop['the'][1]) # Needed with antenna patterns.
#   print('Elev angles: %.3g, %.3g deg' % (el1,el2))
    # Retrieve lrProfile() results to print.
    dh = itm.prop['dh']
    dist_m = itm.prop['dist']
    dla = itm.propa['dla']
    dlsa = itm.propa['dlsa']
    dx = itm.propa['dx']
    ens = itm.prop['ens']
    he = itm.prop['he']

    fs_dB = itm.freespace_dB(1e6*f_MHz, dist_m)
    print('%21s%8.1f%8.1f M' % ('EFFECTIVE HEIGHTS', he[0], he[1]))
    print('%21s%12.0f M' % ('TERRAIN, DELTA H', dh))
    print('')
    gammaE = itm.prop['gme'] # Effective earth curvature within 1 km of surface.
    q = gammaA/gammaE
    print('   POL=%d, EPS=%3.0f, SGM=%6.3f S/M' % (ipol, eps, sgm))
    print('   CLIM=%d, NS=%4.0f, K=%6.3f' % (klim, ens, q))
    print('   PROFILE- NP=%4d, XI=%6.3f KM' % (npt, xkm))
    print('')

    q = dist_m - dla
    q = max(0, q-0.5*pfl[1]) - max(0, -q-0.5*pfl[1])
    if q < 0:
        print('      A LINE-OF-SIGHT PATH')
    else:
        if q == 0:
            print('      A SINGLE HORIZON PATH')
        else: # q > 0
            print('      A DOUBLE-HORIZON PATH')
        if dist_m <= dlsa:
            print('      DIFFRACTION IS THE DOMINANT MODE')
        elif dist_m <= dx:
            print('      TROPOSPHERIC SCATTER IS THE DOMINANT MODE')

#
#           COMPUTE AND PRINT VALUES
    print('')
    print('   ESTIMATED QUANTILES OF BASIC TRANSMISSION LOSS (DB)')
    print('      FREE SPACE VALUE-%7.1f DB' % fs_dB)
    print('')
    print('         RELIA-    WITH CONFIDENCE')
    print('         BILITY ', end='')
    for jc in range(nc):
        print('%8.1f' % qc[jc], end='')
    print('\n')
    xlb = np.zeros(nc)
    for jr in range(nr):
        print('    %10.1f' % qr[jr], end='')
        for jc in range(nc):
            xlb[jc], _ = fs_dB + itm.aVar(zr[jr], 0, zc[jc])
#            print('aVar(%.1f, %.1f, %.1f) ' % (zr[jr],0,zc[jc]), end='')
            if jc == 0:
                print('%10.1f' % xlb[jc], end='')
            else:
                print('%8.1f' % xlb[jc], end='')
        print('')

#
#           PRINT ERROR MESSAGES
    kwx = itm.propv['err']
    if kwx == 0: # No error.
        pass
    elif kwx == 1:
        print('\n')
        print('   **WARNING- SOME PARAMETERS ARE NEARLY OUT OF RANGE.')
        print('     RESULTS SHOULD BE USED WITH CAUTION.')
    elif kwx == 2:
        print('\n')
        print('   **NOTE- DEFAULT PARAMETERS HAVE BEEN SUBSTITUTED')
        print('     FOR IMPOSSIBLE ONES.')
    elif kwx == 3:
        print('\n')
        print('   **WARNING- A COMBINATION OF PARAMETERS IS OUT OF RANGE.')
        print('     RESULTS ARE PROBABLY INVALID.')
    elif kwx == 4:
        print('\n')
        print('   **WARNING- SOME PARAMETERS ARE OUT OF RANGE.')
        print('     RESULTS ARE PROBABLY INVALID.')

    if wqit:
        break # End program.
