#!/usr/bin/env python

#      This code is taken from "A Guide to the Use of the ITS Irregular
#      Terrain Model in the Area Prediction Mode," by G.A. Hufford,
#      A.G. Longley, and W.A. Kissick.  Run as follows:
#
#      demoArea < demoAreaData
#
#      File demoAreaDat follows, where of course leading 'C' is removed,
#      but not leading spaces:
#== = demoAreaDat follows == =
# 1
#TEST PROBLEM 1, QKAREA
# 2        0.        100.      10.       400.      25
# 32       70
# 4        10.       50.       90.       95
#870110    400.      10.       1.        200
# 1
#TEST PROBLEM 2, QKAREA
#861       25
#99
# 1
#TEST PROBLEM 4, QKAREA
# 31       10
# 57       30.       350.      300.      25.       0.02
# 60000    1200.     5.        1000
#88
#
#== = demoAreaDat ends after required blank line == =
#
#      Mike Markowski, mike.ab3ap@gmail.com
#      Dec 2022

#     PROGRAM QKAREA
#        *QUICK AREA*
#        TO ILLUSTRATE THE USE OF THE LONGLEY-RICE MODEL
#          IN THE AREA PREDICTION MODE
#
#        INPUT IS IN 10-COL FIELDS, THE FIRST OF WHICH IS
#          A SEQUENCE OF DIGITS
#        IN PARTICULAR,
#          COL 1 IS THE *EXECUTE* COLUMN--A NON-ZERO DIGIT
#             WILL FORCE OUTPUT
#          COL 2 INDICATES THE CARD TYPE-C~
#
#                          COL
#                          12      11,..
#            STOP-         XO         (OR A BLANK CARD)
#            TITLE-        X1         (NEXT CARD HAS 60-COL TITLE)
#            DISTANCES-    X2      D0,D1,DS1,D2,DS2
#            RELIABILITY-  X3V     QT,QL
#            CONFIDENCE-   X4      QC1,QC2,..
#            ENVIRONMENT-  X5C     DH,NO,ZS,EPS,SGM
#            SYSTEM-       X6NPSS  FMHZ,HG1,HG2
#            (ALTERNATE)   X7NPSS FMHZ, HG 1, HG2, DH, NS, EPS, SGM
#            EXECUTE-      X8
#            RESET-        X9
#

from scipy.stats import norm
import itm
import numpy as np
import sys

itm.dictionaries() # Create dictionaries.

qc = np.zeros(7)
zc = np.zeros(7)
xin = np.zeros(7)

gma = 157e-9
db = 8.685890
akm = 1000
wcon = True
fmhz = 100
hg = np.array([3, 3])
dh = 90
en0 = 301
zsys = 0
eps = 15
sgm = 0.005
ipol = 1
kst = ['r', 'r']
klim = 5-1
itm.mdvarSet(itm.mBcast)
nc = 3
qc[:3] = np.array([50, 90, 10])
zc[:3] = norm.isf(qc[:3]/100) # 0., -1.28155, 1.28155
qt = ql = 50
zt = zl = 0
d0 = 10
ds = 10
dsc = 50
nd = 22
ndc = 15
wtl = False

while True:
    # Replicate Fortran 6I1 format
    jin = np.zeros(6, dtype=int)
    try:
        line = input()
    except EOFError:
        sys.exit(0)
    for i in range(min(6, len(line))):
        jin[i] = int(line[i]) if line[i].isnumeric() else 0
    line = line[10:] # Chop off 6I1,4X
    wcon = (jin[0] == 0) # Is 1st column zero
    jq = jin[1] # 2nd column of input line

    # Replicate Fortran 7F10.0 format.
    n = int(np.ceil(len(line)/10)) # Ten or fewer F10.0 floats follow.
    if jq > 0:
        xin = np.zeros(7)
        for i in range(n):
            xin[i] = float(line[:10].replace(' ', '0'))
            line = line[10:]

    if jq == 0:
         sys.exit(0)

    elif jq == 1:
        itl = input()
        wtl = True

    elif jq == 2:
        xin[0] = max(0, xin[0])
        q = xin[1] - xin[0]
        if q <= 0:
            if xin[0] != 0:
                d0 = xin[0]
                ds = 0
                dsc = 0
                no = 1
                ndc = 0
        else:
            # xin [start, stop0, step0, stop1, step1]
            if xin[2] <= 0:
                xin[2] = max(1, int(q/20 + 0.5)) # 20 steps if unspecified.
            if xin[0] <= 0:
                xin[0] = xin[2] # 1st distance is step size if unspecified.
            d0 = xin[0] # Starting distance.
            ds = xin[2] # Distance step.
            dsc = ds
            nd = int(max(0, xin[1] - xin[0])/ds + 1.75)
            ndc = 0
            if xin[3] > xin[1]:
                if xin[4] <= 0:
                    xin[4] = 5*xin[2]
                dsc = xin[4]
                jq = (xin[3] - xin[1])/dsc + 0.75
                ndc = nd
                nd += jq

    elif jq == 3: # Set QT, QL.
        # Set mode of variability: 0 (single msg), 1 (accidental),
        # 2 (mobile), 3 (broadcast).
        match jin[2]:
            case 0: itm.mdvarSet(itm.mSingle)
            case 1: itm.mdvarSet(itm.mRandom)
            case 2: itm.mdvarSet(itm.mMobile)
            case 3: itm.mdvarSet(itm.mBcast)
            case _: itm.mdvarSet(itm.mBcast)
        qt = 50 # Default quantile of time.
        ql = 50 # Default quantile of location.
        zt = 0
        zl = 0
        if xin[0] > 0:
            qt = xin[0]
            zt = norm.isf(qt/100)
        elif xin[1] > 0:
            ql = xin[1]
            zl = norm.isf(ql/100)

    elif jq == 4: # Set confidence QC1, CQ2.
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
        if jin[2] > 0:
            itm.propv['klim'] = klim = jin[2]-1
            itm.avarRecalc(itm.lClim) # Tell aVar() that climate changed.
        if xin[0] >= 0:
            itm.prop['dh'] = dh = xin[0]
        if xin[1] > 0:
            en0 = xin[1]
            zsys = xin[2]
        if xin[3] > 0:
            eps = xin[3]
            sgm = xin[4]

    elif jq == 6:
        if jin[2] != 1:
            ipol = min(jin[3], 1)
            siting = min(jin[4], 2)
            kst[0] = ['r', 'c', 'v'][siting]
            siting = min(jin[5], 2)
            kst[1] = ['r', 'c', 'v'][siting]
        if xin[0] > 0:
            fmhz = xin[0]
        if xin[1] > 0:
            itm.prop['hg'][0] = hg[0] = xin[1]
        if xin[2] > 0:
            itm.prop['hg'][1] = hg[1] = xin[2]

    elif jq == 7:
        if jin[2] != 1:
            ipol = min(jin[3], 1)
            siting = min(jin[4], 2)
            kst[0] = ['r', 'c', 'v'][siting]
            siting = min(jin[5], 2)
            kst[1] = ['r', 'c', 'v'][siting]
        if xin[0] > 0:
            fmhz = xin[0]
        if xin[1] > 0:
            itm.prop['hg'][0] = hg[0] = xin[1]
        if xin[2] > 0:
            itm.prop['hg'][1] = hg[1] = xin[2]
        if xin[3] >= 0:
            itm.prop['dh'] = dh = xin[3]
        if xin[4] > 0:
            en0 = xin[4]
            zsys = 0
        if xin[5] > 0:
            eps = xin[5]
            sgm = xin[6]

    elif jq == 8:
        wcon = False

    elif jq == 9:
        fmhz = 100
        itm.prop['hg'][0] = hg[0] = 3
        itm.prop['hg'][1] = hg[1] = 3
        itm.prop['dh'] = dh = 90
        en0 = 301
        zsys = 0
        eps = 15
        sgm = 0.005
        ipol = 1
        kst[0] = 'r'
        kst[1] = 'r'
        klim = itm.cContTemp
        itm.mdvarSet(itm.mBcast) # Broadcast mode.
        itm.avarRecalc() # Full aVar() recalc.
        nc = 3
        qc[:3] = np.array([50, 90, 10])
        zc[:3] = norm.isf(qc[:3]/100) # 0., -1.28155, 1.28155
        qt = ql = 50
        zt = zl = 0
        d0 = 10
        ds = 10
        dsc = 50
        nd = 22
        ndc = 15
        wtl = False

    if wcon:
        continue

#
#       EXECUTION
#
    itm.errorSet(itm.eNone) # No error.
    pol = 'h' if ipol==0 else 'v' # Antenna polarization.

    itm.propv['klim'] = klim
    itm.lrPrep(fmhz, zsys, en0, pol, eps, sgm)
    ens = itm.prop['ens'] # surface refractivity, N-units.
    gme = itm.prop['gme'] # effective earth curvature.

    itm.lrArea(kst)
    he = itm.prop['he']   # Antenna effective heights.
#
#     WRITE HEADING
#
    if not wtl:
        print('   AREA PREDICTIONS FROM THE LONGLEY-RICE MODEL, VERSION 1.2.2')
    else:
        print('   %s' % itl)
    print('')
    print('')
    print('%21s%12.0f MHZ' % ('FREQUENCY', fmhz))
    print('%21s%8.1f%8.1f M' % ('ANTENNA HEIGHTS', hg[0], hg[1]))
    print('%21s%8.1f%8.1f M  (SITING=%c,%c)'
        % ('EFFECTIVE HEIGHTS', he[0], he[1], kst[0].upper(), kst[1].upper()))
    print('%21s%12.0f M' % ('TERRAIN, DELTA H', dh))
    print('')

    q = gma/gme
    print('   POL=%d, EPS=%3.0f, SGM=%6.3f S/M' % (ipol, eps, sgm))
    print('   CLIM=%d, N0=%4.0f, NS=%4.0f, K=%6.3f'
        % (klim, en0, ens, q))
    print('')

    if itm.mdvarGet(itm.mSingle):
        print('   SINGLE-MESSAGE SERVICE')
    elif itm.mdvarGet(itm.mRandom):
        print('   ACCIDENTAL SERVICE')
        print('        %5.1f PER CENT TIME AVAILABILITY' % qt)
    elif itm.mdvarGet(itm.mMobile):
        print('   MOBILE SERVICE')
        print('        REQUIRED RELIABILITY-%5.1f PER CENT' % qt)
    elif itm.mdvarGet(itm.mBcast):
        print('   BROADCAST SERVICE')
        print('        REQUIRED RELIABILITY-%5.1f PER CENT TIME' % qt)
        print('%29s%5.1f PER CENT LOCATIONS' % (' ', ql))
    print('')

#
#       COMPUTE AND PRINT VALUES
    print('   ESTIMATED QUANTILES OF BASIC TRANSMISSION LOSS(DB)')
    print('')
    print('       DIST     FREE    WITH CONFIDENCE')
    print('        KM     SPACE', end='')
    for jc in range(nc):
        print('%8.1f' % qc[jc], end='')
    print('\n')
    dt = ds
    d = d0
    for jd in range(1,int(nd)+1):
        itm.avarRecalc(itm.lDist)
        itm.lrProp(d*akm) # Calculate Aref, reference attenuation.
        err = itm.propv['err']
        fs = itm.freespace_dB(1e6*fmhz, d*akm)
        xlb = np.zeros(nc)
        for jc in range(nc):
            xlb[jc], _ = fs + itm.aVar(zt, zl, zc[jc])
        print('  %9.1f%9.1f%9.1f' % (d, fs, xlb[0]), end='')
        for jc in range(1, nc):
            print('%8.1f' % xlb[jc], end='')
        print('')
        if jd == ndc:
            dt = dsc
        d += dt

#        PRINT ERROR MESSAGES
    if itm.errorGet(itm.eRange2):
        print('\n **WARNING- SOME PARAMETERS ARE OUT OF RANGE.')
        print('  RESULTS ARE PROBABLY INVALID.')
    elif itm.errorGet(itm.eRange):
        print('\n  **WARNING- A COMBINATION OF PARAMETERS ', end='')
        print('IS OUT OF RANGE.')
        print('    RESULTS ARE PROBABLY INVALID.')
    elif itm.errorGet(itm.eSubst):
        print('\n   **NOTE- DEFAULT PARAMETERS HAVE BEEN SUBSTITUTED')
        print('    FOR IMPOSSIBLE ONES.')
    elif itm.errorGet(itm.eEdge):
        print('\n   **WARNING- SOME PARAMETERS ARE NEARLY OUT OF RANGE.')
        print('    RESULTS SHOULD BE USED WITH CAUTION.')
