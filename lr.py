#!/usr/bin/env python

# This program is based on code from "A Guide to the Use of the ITS
# Irregular Terrain Model in the Area Prediction Mode," by G.A. Hufford,
# A.G. Longley, and W.A. Kissick.  Run as follows:
#
# lr.py < lrData
#
# Mike Markowski, mike.ab3ap@gmail.com
# May 2023

from scipy.stats import norm
import itm
import numpy as np
import sys

def config():
    '''Provide an configuration dictionary with default values.  Some are
    used only in area mode, some only in profile mode, but all are created.
    '''

    d = {}

    # Climate indices:
    # cEquatorial    Equatorial.
    # cContSubtrop   Continental subtropical.
    # cMariSubtrop   Marine subtropical.
    # cDesert        Desert.
    # cContTemp      Continental temperate.
    # cMariTempLand  Maritime temperate, overland.
    # cMariTempSea   Maritime temperate, oversea.
    d['climate'] = itm.cContTemp    # Continental temperate.
    d['qSit'] = np.array([50, 90, 10]) # Confidences, quantiles of situations.
    d['deltaH'] = 90    # float, terrain irregularity parameter.
    d['dist'] = []      # float, distance from a station.
    d['epsilon'] = 15   # float, relative ground permittivity, F/m.
    d['freq'] = 100     # float, frequency of wave, MHz.
    d['hgtRx'] = 3      # float, height of antenna structure above surface, m.
    d['hgtTx'] = 3      # float, height of antenna structure above surface, m.
    # Transmission mode, aka mode of variability:
    # mSingle  Single message mode.
    # mRandom  Individual/Random receiver of broadcast.
    # mMobile  Mobile mode.
    # mBcast   Broadcast mode.
    d['modeVar'] = itm.mBcast # Broadcast mode.
    d['N0'] = 301       # Sea level refractivity.
    d['Ns'] = 0         # Surface level refractivity.
    d['pathLen'] = 0    # Sea level distance between 2 stations.
    d['pfl'] = []       # Terrain profile, m.
    d['pol'] = 'v'      # Vertical.
    d['pTitle'] = ''    # Profile title.
    d['qTime'] = np.array([50]) # Quantile of time.
    d['qLoc'] = np.array([50])  # Quantile of location.
    d['qSit'] = np.array([50])  # Quantile of time.
    d['sigma'] = 0.005  # float, relative ground conductivity, S/m.
    d['sitingRx'] = 'r' # char, stn siting r(andom), c(areful), v(ery careful).
    d['sitingTx'] = 'r' # char, stn siting r(andom), c(areful), v(ery careful).
    d['title'] = ''     # Title of test or experiment.
    d['zSys'] = 0       # System elevation above sea level, m.
    return d

def errorMsg():
    e, w = itm.errors()
    if len(w) > 0:
        print('\n   **Warning:')
        for i in range(len(w)):
            print('     %s' % w[i])
        print('     Results should be used with caution.')
    if len(e) > 0:
        print('\n   **Error:')
        for i in range(len(e)):
            print('     %s' % e[i])
        print('    Results are probably invalid.')

def hdrArea(d):
    ens = itm.prop['ens'] # surface refractivity, N-units.
    gme = itm.prop['gme'] # effective earth curvature.
    he  = itm.prop['he']  # Antenna effective heights.

    if d['title'] == '':
        print('   Area Predictions from the Longley-Rice Model, '
            'version %s' % itm.version)
    else:
        print('   %s' % d['title'])
    print('')
    print('%21s%12.3f MHz' % ('Frequency', d['freq']))
    print('%21s%8.1f%8.1f m' % ('Antenna heights', d['hgtTx'], d['hgtRx']))
    print('%21s%8.1f%8.1f m  (siting=%c,%c)' % ('Effective heights',
        he[0], he[1], d['sitingTx'], d['sitingRx']))
    print('%21s%12.0f m' % ('Terrain, delta H', d['deltaH']))
    print('')

    K = 1/(itm.rEarth_m*gme)
    print('   Polarity=%c, epsilon=%.0f, sigma=%.3f S/m'
        % (d['pol'], d['epsilon'], d['sigma']))
    print('   Climate=%d, N0=%.0f, Ns=%.0f, K=%.3f'
        % (d['climate'], d['N0'], ens, K))
    print('')

    if d['modeVar'] == itm.mSingle:
        print('   First-try success service')
    elif d['modeVar'] == itm.mRandom:
        print('   Random receiver service')
        print('        %.1f percent time availability' % (100*d['qTime'][0]))
    elif d['modeVar'] == itm.mMobile:
        print('   Mobile service')
        print('        Required reliability %.1f percent' % (100*d['qTime'][0]))
    elif d['modeVar'] == itm.mBcast:
        print('   Broadcast service')
        print('        Required reliability %.1f percent time'
            % (100*d['qTime'][0]))
        print('%28s%5.1f percent locations' % (' ', (100*d['qLoc'][0])))
    print('')

def hdrPfl(d):
    print('\n')
    if d['title'] != '': # Test title from input.
        print('   %s' % d['title'])
    else: # Default test data.
        print('Link Predictions from the Longley-Rice Model, Version 1.2.2')
    print('')
    if d['pTitle'] != '': # Profile title.
        print('   %s' % d['pTitle'])
    print('')
    print('%21s%12.1f km' % ('Distance', d['pathLen']))
    print('%21s%12.3f MHz' % ('Frequency', d['freq']))
    print('%21s%8.1f%8.1f m' % ('Antenna heights', d['hgtTx'], d['hgtRx']))

def runArea(d):
    '''
    runArea() runs an ITM simulation in area mode.  It naively assumes
    that all values needed exist in the passes dictionary.  If they
    do not, the program will die ungracefully.  The real goal of the
    subroutine is to collect everything needed so that a call can be
    made to lrArea().  The results of that call are modified by calls
    to aVar() and printed in a manner similar to the old Fortran program
    qkarea from ITM v1.2.1.

    Input:
      d (dictionary): dictionary of values read in from an rfl config file.
    Output:
      none, all results to stdout.
    '''

    itm.errorSet(itm.eNone) # No error.
    itm.prop['dh'] = d['deltaH']
    itm.prop['hg'][0] = d['hgtTx']
    itm.prop['hg'][1] = d['hgtRx']
    itm.propv['klim'] = d['climate']
    itm.propv['mdvar'] = d['modeVar']

    conf    = d['qSit']
    epsilon = d['epsilon']
    f_MHz   = d['freq']
    N0      = d['N0']
    pol     = d['pol']
    sigma   = d['sigma']
    zSys    = d['zSys']
    zS      = norm.isf(d['qSit'])
    zL      = norm.isf(d['qLoc'][0])
    zT      = norm.isf(d['qTime'][0])

    itm.lrPrep(f_MHz, zSys, N0, pol, epsilon, sigma)
    itm.lrArea([d['sitingTx'], d['sitingRx']])

    hdrArea(d)
    # Compute and print values.
    print('   Estimated Quantiles of Basic Transmission Loss (dB)\n')
    print('       Dist     Free    With Confidence')
    print('        km     Space ', end='')
    for i in range(len(conf)):
        print('%8.1f' % (100*conf[i]), end='')
    print('')
    print('       ----    -----    ', end='')
    for i in range(len(conf)):
        print('-----   ', end='')
    print('')

    for distRange in d['dist']:
        start, stop, delta = distRange
        for dist in np.arange(start, stop+1, delta):
            itm.lrProp(dist*1e3) # Calculate Aref, reference attenuation.
            fs = itm.freespace_dB(f_MHz*1e6, dist*1e3)
            print('  %9.1f%9.1f ' % (dist, fs), end='')
            for i in range(len(conf)):
                xlb, _ = fs + itm.aVar(zT, zL, zS[i])
                print('%8.1f' % xlb, end='')
            print('')
    errorMsg()
    print('')

def runPfl(d):
    '''
    runPfl() runs an ITM simulation in the profile mode.  It naively
    assumes that all values needed exist in the passes dictionary.  If they
    do not, the program will die ungracefully.  The real goal of the
    subroutine is to collect everything needed so that a call can be made to
    lrProfile().  The results of that call are modified by calls
    to aVar() and printed in a manner similar to the old Fortran program
    the old Fortran program qkpfl from ITM v1.2.2.

    Input:
      d (dictionary): dictionary of values read in from an rfl config file.
    Output:
      none, all results to stdout.
    '''

    hdrPfl(d)
    itm.mdvarSet(itm.mRandom) # Random receiver (of a broadcast) mode.
    itm.mdvarSet(itm.mPfl)    # Point to point mode.

    itm.errorSet(itm.eNone) # No error.
    itm.prop['hg'][0] = d['hgtTx']
    itm.prop['hg'][1] = d['hgtRx']
    itm.propv['klim'] = d['climate']

    conf    = d['qSit']
    epsilon = d['epsilon']
    f_MHz   = d['freq']
    N0      = d['N0']
    Ns      = d['Ns']
    pfl     = np.array(d['pfl'])
    pol     = d['pol']
    rel     = d['qTime']
    sigma   = d['sigma']
    zSys    = d['zSys']
    zS      = norm.isf(d['qSit'])
    zT      = norm.isf(d['qTime'])

    xkm = d['pathLen'] # Length of entire path.
    dkm = xkm/(pfl.size - 1)
    rng = np.arange(pfl.size)*dkm*1e3

    Nsp = Ns
    if Ns <= 0: # Calculate zSys.
        i0 = int(0.1*pfl.size)
        i1 = int(0.9*pfl.size)
        zSys = d['zSys'] = pfl[i0:i1+1].mean()
        Nsp = N0
    itm.lrPrep(f_MHz, zSys, Nsp, pol, epsilon, sigma)
    itm.lrProfile(rng, pfl)
    dh = itm.prop['dh']
    dist_m = itm.prop['dist']
    dla = itm.propa['dla']
    dlsa = itm.propa['dlsa']
    dx = itm.propa['dx']
    Ns = itm.prop['ens']
    he = itm.prop['he']

    fs_dB = itm.freespace_dB(1e6*f_MHz, dist_m)
    print('%21s%8.1f%8.1f m' % ('Effective heights', he[0], he[1]))
    print('%21s%12.0f m' % ('Terrain, delta h', dh))
    print('')
    gammaE = itm.prop['gme'] # Effective earth curvature within 1 km of surface.
    K = 1/(itm.rEarth_m*gammaE)
    print('   Pol=%c, epsilon=%.0f, sigma=%.3f S/m' % (pol, epsilon, sigma))
    print('   Climate=%d, Ns=%.0f, K=%.3f' % (d['climate'], Ns, K))
    print('   Profile: np=%d, delta d=%.3f km\n' % (pfl.size, dkm))
    pmode = ''
    res1, res2 = itm.propMode(pfl)
    match res1:
        case 'l': pmode = 'A line-of-sight path'
        case 's': pmode = 'A single horizon path'
        case 'd': pmode = 'A double horizon path'
    print('       %s' % pmode)
    match res2:
        case 'd': pmode = 'Diffraction is the dominant mode'
        case 't': pmode = 'Tropospheric scatter is the dominant mode'
    if res2 != '': print('       %s' % pmode)

    # Compute and print values.
    print('')
    print('   Estimated quantiles of basic transmission loss (dB)')
    print('      Free space value: %.1f dB' % fs_dB)
    print('')
    print('         Relia-    With Confidence')
    print('         bility', end='')
    for i in range(len(conf)):
        print('%8.1f' % (100*conf[i]), end='')
    print('')
    print('         ------    ', end='')
    for i in range(len(conf)):
        print('-----   ', end='')
    print('')
    for r in range(len(zT)): # Reliability index.
        print('    %10.1f' % (100*rel[r]), end='')
        for c in range(len(zS)): # Confidence index.
            xlb, _ = fs_dB + itm.aVar(zT[r], 0, zS[c])
            fmt = '%10.1f' if c==0 else '%8.1f' # More spacing for 1st.
            print(fmt % xlb, end='')
        print('')

def main():
    ''' This subroutine is the main program and is simple.  It reads the
    conifugation file line by line and stores or updates the variables that
    it reads.  The only exception is that if the keyword 'run' is
    encountered, it calls either runPfl() and runArea() depending which ITM
    mode is selected by the user.  When there are no more lines to be read
    from stdin the program exits.
    '''

    itm.dictionaries() # Create dictionaries.
    d = config()
    newRun = True
    while True:
        try:
            line = input()
        except EOFError:
            sys.exit(0) # Exit after reading entire input file.

        excl = line.find('#')
        if excl != -1: # Ignore comments.
            line = line[:excl].strip()
        line = line.strip()
        if line == '': # Ignore blank lines.
            continue
        field = line.split()
        if field[0] in ['climate']:
            d[field[0]] = int(field[1])
        elif field[0] in ['modeVar']:
            # Convert old style index to new style.
            match int(field[1]):
                case 0: itm.mdvarSet(itm.mSingle)
                case 1: itm.mdvarSet(itm.mRandom)
                case 2: itm.mdvarSet(itm.mMobile)
                case 3: itm.mdvarSet(itm.mBcast)
                case _: print('mdvar error.')
            d['modeVar'] = itm.propv['mdvar']
        elif field[0] in ['deltaH', 'epsilon', 'freq', 'hgtRx', 'hgtTx',
            'N0', 'Ns', 'pathLen', 'sigma', 'zSys']:
            d[field[0]] = float(field[1])
        elif field[0] in ['pol', 'sitingRx', 'sitingTx']:
            d[field[0]] = field[1][0].lower()
        elif field[0] in ['qTime', 'qLoc', 'qSit']:
            d[field[0]] = np.array(list(map(float, field[1:])))/100
        elif field[0] == 'dist':
            if newRun:
                d['dist'] = []
                newRun = False
            d['dist'].append(list(map(float, field[1:])))
        elif field[0] == 'pfl':
            if newRun:
                d['pfl'] = []
                newRun = False
            d['pfl'].extend(list(map(float, field[1:])))
        elif field[0] == 'run':
            newRun = True
            itm.avarRecalc() # Force full itm.aVar() recalculation.
            if field[1] == 'area':
                runArea(d)
            elif field[1] == 'pfl':
                runPfl(d)
            d['title'] = ''
        elif field[0] in ['pTitle', 'title']:
            space = line.find(' ')
            d[field[0]] = line[space+1:].strip()

if __name__ == '__main__':
    main()
