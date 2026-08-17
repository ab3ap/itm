#!/usr/bin/env python

# This program uses itm.py to solve the problem presented in section 7.3 of 
# "A Guide to the Use of the ITS Irregular Terrain Model in the Area
# Prediction Mode," Hufford, Longley, Kissick, NTIA Report 82-100, 1982.
#
# Except for main() at bottom, subroutines are in alphabetical order.
#
# Mike Markowski, mike.ab3ap@gmail.com
# June 2023

from scipy.stats import norm
import itm
import matplotlib.pyplot as plt
import numpy as np

def predict():
    '''Implement problem 7.2 Optimum Television Station Separation, p. 42 in 
    "A Guide to the Use of the ITS Irregular Terrain Model in the Area
    Prediction Mode," Hufford, Longley, Kissick, NTIA Report 82-100, 1982.
    '''

    # Set up values for R3 data set as given in section 7.3 of report.
    itm.dictionaries()       # Prepare ITM.
    f_MHz = itm.prop['f_MHz'] = 410 # UHF tv-like.
    itm.prop['hg'][0] = 275  # m, transmit antenna structure height.
    itm.prop['hg'][1] = 6.6  # m, receive antenna structure height.
    itm.prop['dh'] = 126     # m, roughness factor for hilly terrain.
    itm.prop['klim'] = itm.cContTemp # Continental temperate.
    itm.mdvarSet(itm.mBcast) # Broadcast mode.
    pol = 'h'     # Horizontal antenna polarization.
    epsilon = 15  # F/m, relative ground permittivity.  Average ground.
    sigma = 0.005 # S/m, relative ground conductivity.  Average ground.
    zSys = 1700   # m, system elevation above sea level.
    N0 = 300      # N-units, sea level refractivity (yields Ns=250).
    zT = 0        # qT = 0.5, no time variability.

    itm.lrPrep(f_MHz, zSys, N0, pol, epsilon, sigma)
    itm.lrArea(['r', 'r']) # Recreate figure 9.
#    itm.lrArea(['v', 'v']) # Recreate figure 12.

    plotsX = [] # Start collection of curves.
    plotsY = []
    labels = []
    # All curves exactly match report except A(0.1,0.1).
    for qL,qC in [[0.9,0.9],[0.5,0.9],[0.5,0.5],[0.5,0.1],[0.1,0.1]]:
        zL = norm.isf(qL) # Quantile of location.
        zC = norm.isf(qC) # Quantile of confidence (situation).
        plotX = [] # Start new curve.
        plotY = []
        for dist_km in range(1, 101): # Distance from cliff-top transmitter.
            itm.lrProp(dist_km*1e3) # Longley-Rice propagation calculation.
            aRef_dB,_ = itm.aVar(zT, zL, zC) # Variability of attenuation.
            plotX.append(dist_km)
            plotY.append(aRef_dB)
            if aRef_dB > 40: # Match max atten in Figures 9 and 12.
                break
        plotsX.append(plotX) # Save new curve.
        plotsY.append(plotY)
        labels.append('A(%.1f, %.1f)' % (qL, qC)) # This curve's label.
    return plotsX, plotsY, labels

def present(x, y, lbl):
    # Graph settings to match Figures in report.
    fig = plt.figure(figsize=(8, 4), dpi=120)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title('Predicted Attenuation, R3 Data')
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Attenuation (dB)')
    ax.set_xlim(0, 100)
    ax.set_ylim(-9.9, 40)
    ax.invert_yaxis()
    ax.grid(True)

    for i in range(len(x)):
        plt.plot(x[i], y[i], lw=0.75, label=lbl[i])
    plt.legend()
    plt.show()

def main():
    x, y, lbl = predict()
    present(x, y, lbl)

if __name__ == '__main__':
    main()
