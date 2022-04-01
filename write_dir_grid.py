# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 14:05:53 2022

@author: mc16535
"""

from __future__ import print_function # For progress bar
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.interpolate import griddata
import sys
import time

#%% Progress bar and helper functions

def secsToString(seconds):
    return "{:02d}:{:02d}:{:02d}".format(int(seconds/3600), int((seconds%3600)/60), int(seconds%60))

def formatTo3SF(itn):
    if np.log10(itn) >= 2:
        return "{:4d}".format(int(itn))
    elif np.log10(itn) >= 1:
        return "{:3.1f}".format(itn)
    else:
        return "{:3.2f}".format(itn)

def itnToString(itn):
    if np.log10(itn) >= 18:
        return "{}E".format(formatTo3SF(itn*10**-18))
    elif np.log10(itn) >= 15:
        return "{}P".format(formatTo3SF(itn*10**-15))
    elif np.log10(itn) >= 12:
        return "{}T".format(formatTo3SF(itn*10**-12))
    elif np.log10(itn) >= 9:
        return "{}G".format(formatTo3SF(itn*10**-9))
    elif np.log10(itn) >= 6:
        return "{}M".format(formatTo3SF(itn*10**-6))
    elif np.log10(itn) >= 3:
        return "{}k".format(formatTo3SF(itn*10**-3))
    else:
        return "{:4d} ".format(int(itn))

def updateProgress(startTime, currentItn, totalItn):
    currentTime = time.time()
    timeElapsed = secsToString(currentTime - startTime)
    timeRemains = secsToString((currentTime - startTime) * (float(totalItn) / float(currentItn+1) - 1))
    
    barChars = 50
    percPrint = int(100*float(currentItn+1)/float(totalItn))
    percBarFull = int(barChars*float(currentItn+1)/float(totalItn))
    bar = "".join(["#" for ii in range(percBarFull)]) + "".join([" " for ii in range(barChars-percBarFull)])
    
    print("\r{:3d}%|{}|{}/{} Elap:{} ETA:{}  ".format(percPrint, bar, itnToString(currentItn), itnToString(totalItn), timeElapsed, timeRemains), end='')
    sys.stdout.flush()
    
def line_dir(theta, lamL, lamS, isLong=True):
    """
    Computes line directivity from Miller and Pursey. 

    Parameters
    ----------
    theta : float OR ndarray(float)
        Angles for which line directivity will be calculated.
    lamL : float
        Longitudinal wavelength.
    lamS : float
        Shear wavelength.
    isLong : bool, optional
        Test for whether we're computing longitudinal or shear. The default is
        True.

    Returns
    -------
    directivity : complex OR ndarray(complex)
        Direcitivity values at the desired angles. Has the same shape as theta.

    """
    mu2 = lamL**2 / lamS**2
    sinT2 = np.sin(theta) ** 2
    cosT = np.cos(theta)
    
    if isLong:
        directivity = np.exp(3.0j/4.0 * np.pi) * cosT * (mu2 - 2 * sinT2) / F_0(sinT2, mu2) / np.sqrt(np.pi)
    else:
        directivity = np.exp(5.0j/4.0 * np.pi) * mu2**(5/4) * np.sin(2*theta) * np.sqrt(mu2 * sinT2 - 1, dtype=np.complex128) / F_0(mu2*sinT2, mu2) / np.sqrt(np.pi)
        
    return directivity

def F_0(z2, mu2):
    return (2 * z2 - mu2)**2 - 4 * z2 * np.sqrt(z2 - 1, dtype=np.complex128) * np.sqrt(z2 - mu2, dtype=np.complex128)




#%% Import data

interpN = 500

if len(sys.argv) > 1:
    filename = sys.argv[1]
    if len(sys.argv) == 3:
        try:
            interpN = int(sys.argv[2])
        except ValueError:
            raise ValueError("Invalid value for interpN: {}".format(sys.argv[2]))
else:
    filename = "D_rdm1_45npw.dir_32interp.dat"
    os.chdir(r"C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\Abaqus\FMC Generation\v9\Output\Phase Difference Study\Directivity - read input")
lookupname = filename[:-3]+'lookup.dat'

rmag, pmag = np.zeros((10000)), np.zeros((10000))
dL, dS     = np.zeros((10000), dtype=np.complex128), np.zeros((10000), dtype=np.complex128)
vL, vS     = np.zeros((10000)), np.zeros((10000))

with open(filename, 'r') as file:
    for linecount, line in enumerate(file):
        if linecount == 0:
            # line = line.split(':')
            # freq = np.array(line[1][6:-2].split(' '))
            # freq = freq[freq != '']
            # freq = freq[freq != ' ']
            # freq = freq.astype(float)
            # phasediffs = np.zeros((10000, len(freq)))
            continue
        line = line.replace('inf', '1')
        line = line.split(', ')
        rmag[linecount-1] = float(line[0])
        pmag[linecount-1] = float(line[1])
        dL[linecount-1]   = complex(line[2])
        dS[linecount-1]   = complex(line[3])
        vL[linecount-1]   = float(line[4])
        vS[linecount-1]   = float(line[5])
        # phases = np.array(line[6][1:-2].split(' '))
        # phases = phases[phases != '']
        # phases = phases[phases != ' ']
        # phasediffs[linecount-1, :] = phases.astype(float)
        
pmag = pmag[rmag != 0]
dL   = dL[rmag != 0]
dS   = dS[rmag != 0]
vL   = vL[rmag != 0]
vS   = vS[rmag != 0]
rmag = rmag[rmag != 0]



#%% Interpolate directivities

interpN = 100
p, r = np.linspace(min(pmag), max(pmag), interpN), np.linspace(min(rmag), max(rmag), interpN)
[P, R] = np.meshgrid(p, r)

dL_interp = griddata(
                np.array([pmag, rmag]).T, 
                dL, 
                (P, R), 
                method='cubic',
                fill_value=0.0
            )
dS_interp = griddata(
                np.array([pmag, rmag]).T, 
                dS, 
                (P, R), 
                method='cubic',
                fill_value=0.0
            )

#%% Write out file for use in matlab

R_unravelled, P_unravelled   = R.ravel(), P.ravel()
dL_unravelled, dS_unravelled = dL_interp.ravel(), dS_interp.ravel()
with open(lookupname, 'w') as f:
    f.write("rmag, phimag, dL, dS\n")
    for ii in range(interpN**2):
        f.write("{:.12g}, {:.12g}, {:.12g}, {:.12g}\n".format(R_unravelled[ii], P_unravelled[ii], dL_unravelled[ii], dS_unravelled[ii]))

centre_freq = 5e6
theta = np.linspace(-np.pi/2, np.pi/2, 500)
thetalim = .04
an_dL = line_dir(theta, 6317.012224890781/centre_freq, 3110.2818131859126/centre_freq)
an_dS = line_dir(theta, 6317.012224890781/centre_freq, 3110.2818131859126/centre_freq, isLong=False)
mask = np.abs(theta) < np.pi/6
where_dS_max = np.where(an_dS == np.max(an_dS[mask]))[0][0]
shear_mask = np.logical_or(np.abs(P_unravelled - theta[where_dS_max]) < thetalim, np.abs(P_unravelled + theta[where_dS_max]) < thetalim)
shear_mask = np.logical_and(np.logical_and(shear_mask, R_unravelled > 15e-3), np.abs(dS_unravelled) != 0.0)
max_dS_val = np.mean(np.abs(dS_unravelled[shear_mask]))
max_andS_val = np.max(np.abs(an_dS[mask]))
long_mask = np.abs(P_unravelled) < thetalim
long_mask = np.logical_and(np.logical_and(long_mask, R_unravelled > 15e-3), np.abs(dL_unravelled) != 0.0)
max_dL_val = np.mean(np.abs(dL_unravelled[long_mask]))
max_andL_val = np.max(np.abs(an_dL))
print('{}: meas dL/dS = {:.4f}, an dL/dS = {:.4f}'.format(filename, max_dL_val/max_dS_val, max_andL_val/max_andS_val))

# fig = plt.figure(figsize=())
# plt.scatter(rmag * 10**3, pmag, c=np.abs(dS), s=2, vmin=0, vmax=3.5e-8)
# plt.colorbar()
# plt.show()
# plt.scatter(rmag * 10**3, pmag, c=np.angle(dS), s=2, vmin=-np.pi, vmax=np.pi)
# plt.colorbar()
# plt.show()
# plt.imshow(np.abs(dS_interp.T), extent=[np.min(R) * 10**3, np.max(R) * 10**3, np.min(P), np.max(P)], aspect='auto', origin='upper', vmin=0, vmax=3.5e-8)
# plt.colorbar()
# plt.show()
# plt.imshow(np.angle(dS_interp.T), extent=[np.min(R) * 10**3, np.max(R) * 10**3, np.min(P), np.max(P)], aspect='auto', origin='upper', vmin=-np.pi, vmax=np.pi)
# plt.colorbar()
# plt.show()