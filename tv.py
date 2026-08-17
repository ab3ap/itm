#!/usr/bin/env python

# This program uses itm.py to solve the problem presented in section 7.2 of 
# "A Guide to the Use of the ITS Irregular Terrain Model in the Area
# Prediction Mode," Hufford, Longley, Kissick, NTIA Report 82-100, 1982.
# Equation numbers in comments below reference those in the report.
#
# Except for main() at bottom, subroutines are in alphabetical order.
#
# Mike Markowski, mike.ab3ap@gmail.com
# May 2023

from scipy.stats import norm
import itm
import matplotlib.pyplot as plt
import numpy as np

memQL = {} # Memory of results for deviation().
memQT = {} # Memory of results for deviation().

def coverage(res, Rn, s):
    '''Linearly interpolate results to find percent coverage.

    res[] is a 2 row list of lists.  Each row is [x, Rx], where x is distance
    from desired antenna, and Rx is the D/U ratio at x.  The input res looks
    like

    [[x0, Rx0],
     [x1, Rx1]]

    where Rn is in interval [Rx0,Rx1] (or [Rx1,Rx0] if Rx1 < Rx0).
    Rn's fraction f in [Rx0,Rx1] is used to interpolate a value
    for r in [x0,x1].

    Inputs:
      res (float[2][2]): two lines of dist and D/U results.
      Rn (float): R0 or R1, D/U value to interpolate between res[:,1] values.
      s (float): km, separation between broadcast antennas.
    Outputs:
      rho (float): percent interference-free coverage.
      r (float): km, service range of interference coverage.
    '''
    span = res[1] - res[0]       # 1x2 vector, row1 - row0.
    f = (Rn - res[0][1])/span[1] # Rn fraction between Rx's (D/U's).
    res = res[0] + f*span        # Interpolate row by f, now Rn==res[1].
    # Solve Eq 11 and 17 D/U ratio for x=r, interference-free range.
    r = res[0]                   # Dist r when R(r) = Rn = reqd D/U = res[1].

    rho = 100*2*np.pi/np.sqrt(3)*(r/s)**2 # Eq 12, percent coverage.
    return rho, r                # Coverage, interpolated dist.

def deviation(dist_km, qT, qL):
    '''After itm.lrProp() has been called, this subroutine is called with
    2 of the 3 quantiles of attenuation so that deviations from medians are
    calculated and returned.  Eq 19 is implemented to generate quantiles
    of individual attenuation.  Service range is based on time and location
    and is why this subroutine deals with them.  This subroutine is the
    foundation of the program.

    The definition of quantiles of observations are reversed here.  Quantiles
    of observations are values exceeded for given fractions of time, locations
    and situation.  See discussion under Eq 3 in 1982 report.  If parameters
    qT and qL from the calling code represent (1-qT) and (1-qL), as in the
    case of undesired signals, probabilities are that attenuation does not
    exceed aVar(qT,qL,qS) for >= 1-qL locations for >= 1-qT of the time.

    Inputs:
      dist_km (float): distance from antenna to calculate path loss.
      qT, qL (floats): quantiles, 0 to 1, of time and location.
    Outputs:
      Ao (float): overall median of attenuation, i.e., aVar(0,0,0).
      YT (float): time deviation from median.
      YL (float): location deviation from median.
      pl_dB (float): path loss with confidence qT.
    '''

    # Save isf() results for re-use and significant speed up.
    if qT not in memQT.keys():
        memQT[qT] = norm.isf(qT)
    if qL not in memQL.keys():
        memQL[qL] = norm.isf(qL)
    zT = memQT[qT]
    zL = memQL[qL]
    zS = 0 # isf(qS=0.5) == 0, qS=0.5 always for interference probs.

    # Make ITM recalculate distance based vals.
    itm.lrProp(dist_km*1e3) # Longley-Rice propagation calculation.

    # Find time deviation from median for given location.
    # At least qL locations & qT times, atten does not exceed aVar(zT,zL,zS).
    YT1,_ = itm.aVar(zT, zL, zS) # Choose T after fixing L and S.
    YT0,_ = itm.aVar(0., zL, zS) # Time median.

    # Find location deviation from median.
    # For at least qL of locations, attenuation does not exceed aVar().
    YL1,_ = itm.aVar(0., zL, zS) # Choose L after fixing S.
    YL0,_ = itm.aVar(0., 0., zS) # Location median.

    YT = YT1 - YT0 # Deviation from time median.
    YL = YL1 - YL0 # Deviation from location median.

    Ao,_ = itm.aVar(0., 0., 0.)  # dB, overall median of attenuation.
    itm.mdvarUnset(itm.mIfere)   # Otherwise, aVar sets zS to 0.
    aConf,_ = itm.aVar(0, 0, zT) # Use time quantile as confidence.
    itm.mdvarSet(itm.mIfere)     # Restore interference mode.
#   pl_dB = aConf + 20*np.log10(2*itm.prop['wn']*1e3*dist_km)
    pl_dB = aConf + itm.freespace_dB(1e6*itm.prop['f_MHz'], 1e3*dist_km)
    return Ao, YT, YL, pl_dB

def du(distD, distU, qT, qL):
    '''Find the D/U ratio for a given distance and quantiles of attenuation.

    Inputs:
      distD (float): km, distance receiver is from desired antenna.
      distU (float): km, distance receiver is from undesired antenna.
      qT (float): 0 to 1, quantile of time.
      qL (float): 0 to 1, quantile of location.
    Output:
      Rx (float): dB, desired/undesired ratio at distance distD km.
    '''

    # Find Rxo = D/U for closer offset stations.
    # Quantiles of attenuation for time and location.
    AoD, YTD, YLD, plD = deviation(distD,   qT,   qL) # D desired.
    AoU, YTU, YLU, plU = deviation(distU, 1-qT, 1-qL) # U undesired.
    # Eq 18, pseudo-convolutions, quantiles of time,loc for Rx.
    YTR = np.sign(0.5-qT) * np.sqrt(YTD**2 + YTU**2) # Time.
    YLR = np.sign(0.5-qL) * np.sqrt(YLD**2 + YLU**2) # Location.
    # Eq 17, Rx is D/U desired-to-undesired ratio.
    Rx = 20*np.log10(distU/distD) # Freespace D/U.
    Rx += -AoD + AoU # Overall medians of attenuation.
    Rx += YLR + YTR  # Time/loc quantiles of atten for ratio R0.
    return Rx, plD, plU   # D/U ratio at distD.

def present(tbl):
    '''Present results to user in two ways: text based table and graph.  The
    table uses a '>' to indicate the optimal separation and also uses an '*'
    to show when the farther non-offset antenna determines the range.

    Inputs:
      tbl (float[][5]): table of results to be displayed.  Each row of table
        contains [sep_km, cov_pct, rng_km, du_dB, plD_dB, plU_dB]:
          sep_km (float): broadcast antenna separation in triangular grid.
          cov_pct (float): percent coverage at Grade A service.
          rng_km (float): range that meets minimum D/U.
          plD_dB (float): desired sig path loss with confidence qT.
          plU_dB (float): undesired sig path loss with confidence qT.
        Each row's results are based on the given separation distance and
        du_dB is always either R0 or R1, the required D/U ratios for offset
        and non-offset frequencies.  The goal is to find the optimal separation
        to achieve maximum coverage.
    Output:
      none
    '''

    best = np.where(tbl[:,1]==max(tbl[:,1]))[0][0] # Index of max coverage.
    # Print table of results.
    R1 = np.max(tbl[:,5]) # Max D/U considered.
    print(' %6s %6s %7s %7s %7s %7s'
        % ('s km', 'cov %', 'r km', 'plD dB', 'plU dB', 'D/U dB'))
    for rnum, r in enumerate(tbl):
        print('%s' % '>' if rnum==best else ' ', end='') # Optimal s line.
        print('%6.1f %6.1f %7.1f %7.1f %7.1f %7.1f'
            % (r[0],r[1],r[2],r[3],r[4],r[5]), end='')
        print('%s' % '+' if r[5]==R1 else '') # Indicator when D/U changes.
    print('\nOptimal separation s = %.1f km.' % tbl[best,0])
    print('Service range      r = %.1f km.' % tbl[best,2])
    print('Percent coverage   %.1f%%.' % tbl[best,1])

    # Calculate EIRP using values from Table 7.
    plD_dB = tbl[best,3]         # For optimal separation s.
    kTB_dBW = -138               # Table 7.
    urban_dB = 19                # Table 7.
    snrMin_dB = 30               # Table 7, required S/N.
    rx_dBi = 2                   # Table 7.
    losses_dB = 2                # Table 7.
    multiPath_dB = 6             # Table 6, Rayleigh fading margin.
    uncertainties_dB = 6         # Table 6, margin for uncertainties.
    n_dBW = kTB_dBW + urban_dB   # Thermal end environmental noise.
    sMin_dBW = n_dBW + snrMin_dB # Satisfy minimum S/N.
    s_dBW = sMin_dBW + plD_dB - rx_dBi
    s_dBW += losses_dB + multiPath_dB #+ uncertainties_dB
    print('EIRP required:     %.1f dBW' % s_dBW)

    # Present graph of separation vs coverage.
    fig = plt.figure(figsize=(8, 4), dpi=120)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title('Interference-Free Signal Reception')
    ax.set_xlabel('Station Separation (km)')
    ax.set_ylabel('Area Coverage (%)')
    ax.set_xlim(0, 400) # To match Figure 7 in report.
    ax.set_ylim(0, 30) # To match Figure 7 in report.
    ax.grid(True)
    plt.plot(tbl[:,0], tbl[:,1], lw=0.75)
    plt.show()

def tv():
    '''Implement problem 7.2 Optimum Television Station Separation, p. 42 in 
    "A Guide to the Use of the ITS Irregular Terrain Model in the Area
    Prediction Mode," Hufford, Longley, Kissick, NTIA Report 82-100, 1982.
    '''

    # Set up values as given in report.
    itm.dictionaries()         # Prepare ITM.
    f_MHz = itm.prop['f_MHz'] = 193 # TV channel 10.
    itm.prop['hg'][0] = 305    # m, transmit antenna structure height.
    itm.prop['hg'][1] = 9      # m, receive antenna structure height.
    itm.prop['dh'] = 90        # m, roughness factor for hilly terrain.
    itm.prop['klim'] = itm.cContTemp # Continental temperate.
    itm.prop['ens'] = 301      # N-units, continental temperate refractivity.
    itm.mdvarSet(itm.mBcast) # Broadcast mode.
    itm.mdvarSet(itm.mIfere) # Interference.
    qT = 0.9      # Quantile of time, Grade A tv service problem.
    qL = 0.7      # Quantile of location, Grade A tv service problem.
    pol = 'h'     # Horizontal antenna polarization.
    epsilon = 15  # F/m, relative ground permittivity.
    sigma = 0.005 # S/m, relative ground conductivity.
    zSys = 0      # m, system elevation above sea level (0 sets Ns to N0).
    N0 = 301      # N-units, sea level refractivity.
    # Required desired/undesired ratios:
    R0 = 28 # Offset freqs: 28 dB
    R1 = 45 # Non-offset f: 45 dB.
    itm.lrPrep(f_MHz, zSys, N0, pol, epsilon, sigma)
    itm.lrArea(['r', 'r']) # Siting not given in report.

    # Rx=D/U for closest ring of offset freq antennas.
    tblO = [] # Offset frequency results.
    res = [[],[]] # Save only 2 most recent results.
    i = 0
    for s in range(40, 401, 5): # km, broadcast antenna separation.
        for x in range(1, s): # km, rcvr dist from desired tx antenna.
            i ^= 1 # Alternate 0 and 1 index.
            Rx, plD, plU = du(x, s-x, qT, qL) # D/U ratio.
            res[i] = np.array([x, Rx]) # New results.
            if Rx < R0: # Stop when Rx=D/U falls below R0.
                rho, r = coverage(res, R0, s) # Solve Eq 11 and 12.
                tblO.append([s, rho, r, plD, plU, R0]) # Build table of results.
                break # Found R0 for this s, stop looking.
    tblO = np.array(tblO)

    # Rx=D/U for farther ring of non-offset freq antennas.
    tblN = [] # Non-offset frequency results.
    for s in range(40, 401, 5): # km, broadcast antenna separation.
        for x in range(1, 2*s): # km, rcvr dist from desired tx antenna.
            i ^= 1 # Alternate 0 and 1 index.
            Rx, plD, plU = du(x, np.sqrt(3)*s-x, qT, qL) # D/U ratio.
            res[i] = np.array([x, Rx]) # New results.
            if Rx < R1: # Stop when Rx=D/U falls below R1.
                rho, r = coverage(res, R1, s) # Solve Eq 11 and 12.
                tblN.append([s, rho, r, plD, plU, R1]) # Build table of results.
                break # Found R1 for this s, stop looking.
    tblN = np.array(tblN)

    # Save row i from whichever table has shortest r, most interfered signal.
    # Each tbl row is [separation, coverage, range, D/U].
    # XXX Assumes tblO and tblN are same size.  If not, program dies.  Fix me!
    i = tblO[:,2] < tblN[:,2] # Indices for shortest r of R(R0) and R(R1).
    tblN[i] = tblO[i] # tblN now has rows with shortest r's.
    return tblN

def main():
    result = tv()
    present(result)

if __name__ == '__main__':
    main()
