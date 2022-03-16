# coding=utf_8

# Writes .dat files containing directivities for the specified .odb in the current
# working directory. Note that this script must be called from within Abaqus,
# i.e. at the command prompt, use command "abaqus python write_all_direcitivities.py
# <odb_filename>.odb" in the directory containing the files.

# N.B. This code assumes that the yaml module for Python 2.7.3 is installed at
# "C:\Python27\Lib\site-packages" if running on Windows, or "/user/home/mc16535/Admin/
# PythonModules" if running on Linux (BluePebble Cluster). If the user installs the 
# right version of Python and yaml on Windows, this should automatically work; users on
# the cluster or other Linux distros will have to upload/install yaml and manually
# update the path themselves.

# A timetrace is defined as U2 History Output data.

from __future__ import print_function # For progress bar
import numpy as np
import os
import sys
import time

def nanmean(a):
    """
    Equivalent to np.nanmean() from a later version of numpy. For use in Python 2.7
    """
    a = a.flatten()
    a = a[np.invert(np.isnan(a))]
    return np.mean(a)

def read_settings(filename):
    with open(filename, 'r') as f:
        settings = yaml.safe_load(f)
    return settings

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



#%% Compute time of arrival and frequency spectra from input signal
def getSpectrum(u, r, rmag, phi, time0, time1, spec_window, input_duration=1e-6, v_L=6317.012224890781, v_S=3110.2818131859126, interp_factor=8):
    # Get expected time of arrival of longitudinal and shear components
    long_sample_time = rmag / v_L + input_duration/2
    shear_sample_time = rmag / v_S + input_duration/2
		
    # Time-domain window to isolate signal of interest
    l_window = np.ones((u.shape[1], 1))
    l_window[time0 < long_sample_time - 1.05*input_duration/2] = 0
    l_window[time0 > long_sample_time + 1.05*input_duration/2] = 0
    s_window = np.ones((u.shape[1], 1))
    s_window[time0 < shear_sample_time - 1.05*input_duration/2] = 0
    s_window[time0 > shear_sample_time + 1.05*input_duration/2] = 0
    
    # Radial and azimuthal displacement
    r_perp = np.dot(np.array([[0, -1], [1, 0]]), r)
    dr = np.reshape(np.dot(u.T, r)/rmag, (u.shape[1], 1)) * l_window
    dp = np.reshape(np.dot(u.T, r_perp)/rmag, (u.shape[1], 1)) * s_window
    
    # Spectra
    fdr = np.fft.rfft(dr, axis=0)
    fdp = np.fft.rfft(dp, axis=0)
    fdr = fdr[:int(dr.shape[0] / 2)] * spec_window
    fdp = fdp[:int(dp.shape[0] / 2)] * spec_window
    
    # Fourier interpolate to get higher resolution time of propagation
    idr = np.fft.ifft(np.append(interp_factor * fdr, np.zeros((fdr.shape[0]*(interp_factor-1), fdr.shape[1]))), axis=0)
    idp = np.fft.ifft(np.append(interp_factor * fdp, np.zeros((fdp.shape[0]*(interp_factor-1), fdp.shape[1]))), axis=0)
		
    # Get the time of the maximum in the vicinity of the expected time of arrival, to compute actual
    # wave speed.
    time_L = time1[np.where(np.abs(idr) == 
     							np.max(np.abs(idr)))[0][0]] - input_duration/2
    time_S = time1[np.where(np.abs(idp) == 
    							np.max(np.abs(idp)))[0][0]] - input_duration/2
    
    return time_L, time_S, fdr, fdp

def line_dir(theta, lamL, lamS, c44=70e9/(2*(1+0.34)), isLong=True):
    """
    Computes line directivity for a point source from Miller and Pursey.
    Consistent with PW's fn_line_contact_directivity.m function in Matlab, but
    also includes shear modulus as directivity relates force to displacement.

    Parameters
    ----------
    theta : float OR ndarray(float)
        Angles for which line directivity will be calculated.
    lamL : float
        Longitudinal wavelength.
    lamS : float
        Shear wavelength.
    c44 : float, optional
        Shear modulus. Default value assumes E = 70e9 Pa and ν = 0.34.
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
        directivity = np.exp(3.0j/4.0 * np.pi) * np.sqrt(2/np.pi) / c44 * cosT * (mu2 - 2 * sinT2) / F_0(sinT2, mu2)
    else:
        directivity = np.exp(5.0j/4.0 * np.pi) * np.sqrt(2/np.pi) / c44 * mu2**(5/4) * np.sin(2*theta) * np.sqrt(mu2 * sinT2 - 1, dtype=np.complex128) / F_0(mu2*sinT2, mu2)
        
    return directivity

def F_0(z2, mu2):
    return (2 * z2 - mu2)**2 - 4 * z2 * np.sqrt(z2 - 1, dtype=np.complex128) * np.sqrt(z2 - mu2, dtype=np.complex128)



#%% Start. Determine if we are running from command line or not.

try:
    assert len(sys.argv) == 2
    odb_name = sys.argv[1]
    filename = odb_name[:-3]+'dir.dat'
    isCommand = True
    from odbAccess import openOdb
    if "linux" in sys.platform:
        yaml_path = r"/user/home/mc16535/Admin/PythonModules"
        sys.path.append(yaml_path)
    else:
        yaml_path = r"C:\Python27\Lib\site-packages"
        sys.path.append(yaml_path)
# Running from console. yaml must be installed; do not worry about path.
# Save to same directory as write_all_directivities.py, unless running in
# IDE console.
except AssertionError:
    odb_name = 'D_rdm1_45npw.dat'
    filename = odb_name[:-3]+'dir_.dat'
    os.chdir(r"..\AbaqusInputFileGeneration - Output\v9\Output\Phase Difference Study\Directivity - read input")
    isCommand = False
    import matplotlib.pyplot as plt
    
try:
    import yaml
except ModuleNotFoundError:
    raise ModuleNotFoundError('Module "yaml" expected at location {}'.format(yaml_path))

t0 = time.time()

inpname = odb_name[:-3]+'inp'
datname = odb_name[:-3]+'tt.dat'
yamlname = odb_name[:-3]+'yml'

# Input parameters
settings = read_settings(yamlname)

# Input parameters
interp_factor = 16
centre_freq = settings['probe']['freq']
no_cycles = settings['probe']['cycles']
input_duration = no_cycles / centre_freq
density = settings['material']['density']
modulus = settings['material']['modulus']
poisson = settings['material']['poisson']
c44 = modulus / (2. * (1. + poisson))



#%% Read input file to get coordinate locations.

measNodeIdxs = []
measNodeCoords = np.zeros((0, 3))
with open(inpname, 'r') as f:
    # Start by getting the nodes associated with MeasureSet
    isMeasureSet = False
    for line in f:
        # If we're at the start of MeasureSet
        if '*Nset, nset=MeasureSet' in line:
            isMeasureSet = True
            continue
        # If we're at the end of MeasureSet
        if isMeasureSet and ('elset=' in line or 'nset=' in line):
            isMeasureSet = False
            break
        # If we're in the MeasureSet
        if isMeasureSet:
            line = line.split(',')
            measNodeIdxs.append(int(line[0]))

measNodeIdxs = np.unique(measNodeIdxs)
# measNodeIdxs = np.delete(measNodeIdxs, range(1, 8000))

with open(inpname, 'r') as f:
    # Now we have all of the nodes associated with MeasureSet, get their coordinates
    isNodeCoords = False
    for line in f:
        # If we're at the start of the node coordinates
        if '*Node' in line:
            isNodeCoords = True
            continue
        # If we're at te end of the node coordinates
        if measNodeCoords.shape[0] == len(measNodeIdxs):
            isNodeCoords = False
            break
        # If we're in the node coordinates
        if isNodeCoords:
            line = line.split(',')
            if int(line[0]) in measNodeIdxs:
                coords = np.reshape([int(line[0]), float(line[1]), float(line[2])], (1,3))
                measNodeCoords = np.append(measNodeCoords, coords, axis=0)

# Radial distances
r = np.transpose(measNodeCoords[:, 1:])
rmag = np.zeros(r.shape[1])
for ii in range(r.shape[1]):
    rmag[ii] = np.linalg.norm(r[:, ii])

minR = np.where(rmag == min(rmag))[0][0]



#%% Read in timetraces from .odb or .tt.dat, and get input signal as amplitude.

# Read timetraces.
if ".odb" not in odb_name:
    assert ".tt.dat" in datname
    print("Reading from {}".format(datname))
    t1 = time.time()
    with open(datname, 'r') as f:
        kk = 0
        for count, line in enumerate(f):
            line = line.split(',')
            if count == 0:
                timevec = np.array([float(val) for val in line[1:]])
                All_Tts = np.zeros((3, timevec.shape[0], len(measNodeIdxs)))
            elif count == 1:
                input_signal = np.array([float(val) for val in line[1:]])
            elif count%2 == 0:
                if abs(rmag[kk] - float(line[0])) < 1e-8:
                    All_Tts[0, :, kk] = timevec
                    All_Tts[1, :, kk] = np.array([float(val) for val in line[1:]])
                else:
                    continue
            else:
                if abs(np.arcsin(r[0, kk] / rmag[kk]) - float(line[0])) < 1e-8:
                    All_Tts[2, :, kk] = np.array([float(val) for val in line[1:]])
                    kk += 1
                # else:
                #     continue 
                # kk += 1
            if "win" in sys.platform and (kk+1)%int(len(measNodeIdxs)/100) == 0:
                updateProgress(t1, kk+1, len(measNodeIdxs))
            elif "linux" in sys.platform and (kk+1)%int(len(measNodeIdxs)/5) == 0:
                updateProgress(t1, kk+1, len(measNodeIdxs))
            if kk == len(measNodeIdxs):
                break
        print('')
    
    amplitude = np.zeros((timevec.shape[0], 2))
    amplitude[:, 0] = timevec
    amplitude[:, 1] = input_signal
    
    if measNodeIdxs[0] == 1:
        ignoreNode0 = True
    else:
        ignoreNode0 = False

else:
    odb = openOdb(odb_name, readOnly=True)
    
    region = odb.steps['Step-1'].historyRegions
    keys = region.keys()
    
    # Get the input signal.
    if 'CF2' in region[keys[0]].historyOutputs.keys():
        amplitude = np.array(region[keys[0]].historyOutputs['CF2'].data)
        ignoreNode0 = True
        
    else:
        readAmplitude = np.zeros((0, 2))
        with open(inpname, 'r') as f:
            inAmp = False
            inLoad = False
            for line in f:
                if "Amplitude" in line:
                    inAmp = True
                elif inAmp:
                    if "Material" in line:
                        inAmp = False
                        continue
                    line = line.split(',')
                    readAmplitude = np.append(readAmplitude, np.reshape([float(line[0]), float(line[1])], (1,2)), axis=0)
                else:
                    if "Cload" in line:
                        inLoad = True
                    elif inLoad:
                        inLoad = False
                        ampMult = float(line.split(',')[2])
                    else:
                        continue
                    
        readAmplitude[:,1] *= ampMult
        ignoreNode0 = False
        
        time0 = np.array(region[keys[1]].historyOutputs['U1'].data)[:, 0]
        amplitude = np.zeros((time0.shape[0], 2))
        amplitude[:,0] = time0
        amplitude[:,1] = np.interp(time0, readAmplitude[:,0], readAmplitude[:,1], right=0)



#%% Compute input spectrum and get working parameters.
        
time0 = np.array(amplitude[:, 0])
time_pts = amplitude.shape[0]
time_step = time0[1] - time0[0]

fft_pts = time_pts
fstep = 1.0/(time0[-1])
freq = ((np.array(range(int(fft_pts/2)))) * fstep)
sig_window = np.ones((time0.shape[0], 1))
sig_window[time0 > input_duration] = 0
spec_window = np.ones((freq.shape[0], 1))
spec_window[freq > 2*centre_freq] = 0

in_time_sig = np.reshape(amplitude[:, 1], (amplitude.shape[0], 1)) * sig_window
in_freq_spec = np.fft.fft(in_time_sig, axis=0)
in_freq_spec = in_freq_spec[:int(in_time_sig.shape[0] / 2)] * spec_window

# Find out where the spectrum is non-zero for computing directivities later.
nonZeros = np.where((np.abs(in_freq_spec) > max(np.abs(in_freq_spec))/100))[0]
idx1, idx2 = nonZeros[0], nonZeros[-1]

inv_in_freq_spec = np.diag(1.0/np.squeeze(in_freq_spec[idx1:idx2]))

# Fourier interpolate input spectrum to get duration of input signal.
interp_freq_spec = np.zeros((in_freq_spec.shape[0]*interp_factor, 1), dtype=np.complex128)
interp_freq_spec[:in_freq_spec.shape[0]] = interp_factor * in_freq_spec
idu = np.fft.ifft(np.append(interp_factor * in_freq_spec, np.zeros((in_freq_spec.shape[0]*(interp_factor-1), in_freq_spec.shape[1]))), axis=0)
input_duration = 2 * np.linspace(time0[0], time0[-1], idu.shape[0])[np.where(np.abs(idu) == np.max(np.abs(idu)))[0][0]]
        
time1 = np.linspace(time0[0], time0[-1], idu.shape[0])

# Find the time taken for the wave to propagate to the node which is nearest to
# the origin, for calculating wave velocity at each node. Do this for the nearest
# node as error will decrease as radial distance from this node increases, and
# we have fewer nodes nearer to the origin. Using the first node (which has close
# to maximum r) will mean that there will be a larger number of nodes with significant
# error.
if ignoreNode0:
    k1 = minR+1
else:
    k1 = minR

nodeIdx = measNodeIdxs[minR]
if ".odb" in odb_name:
    try:
        assert int(keys[k1].split('.')[1]) == nodeIdx 
    except AssertionError:
        print("Input file node {} does not match odb file node {}".format(keys[k1].split('.')[1], nodeIdx))
        raise ValueError
    u1 = np.array(region[keys[k1]].historyOutputs['U1'].data)
    u2 = np.array(region[keys[k1]].historyOutputs['U2'].data)
else:
    u1 = np.array([All_Tts[0, :, k1], All_Tts[1, :, k1]]).T
    u2 = np.array([All_Tts[0, :, k1], All_Tts[2, :, k1]]).T

u = np.array([u1[:, 1], u2[:, 1]])
phi = np.arcsin(r[0, minR] / rmag[minR])
ref_tL, ref_tS, _, _ = getSpectrum(u, r[:,minR], rmag[minR], phi, time0, time1, spec_window, input_duration=input_duration, interp_factor=interp_factor)
refR = rmag[minR]



#%% Compute directivities

dir_L = np.zeros(measNodeIdxs.shape[0], dtype=np.complex128)
dir_S = np.zeros(measNodeIdxs.shape[0], dtype=np.complex128)
	
t1 = time.time()
print("Time to loop start: {}".format(secsToString(t1 - t0)))

with open(filename, 'w') as f:
    f.write('rmag, phi, C dL, C dS, vL, vS, dL_a - dL_m, dS_a - dS_m')
    f.write('{}\n'.format(np.array_str(freq[idx1:idx2], precision=1, max_line_width=20*(idx2-idx1))))
    for kk in range(measNodeIdxs.shape[0]):
        
        ################################# REMOVE INDENTATIONS WITHIN HERE ##################################
        if ignoreNode0:
            k1 = kk+1
        else:
            k1 = kk
        
        nodeIdx = measNodeIdxs[kk]
        
        # Displacements
        if ".odb" in odb_name:
            try:
                assert int(keys[k1].split('.')[1]) == nodeIdx 
            except AssertionError:
                print("Input file node {} does not match odb file node {}".format(keys[kk].split('.')[1], nodeIdx))
                raise ValueError
            u1 = np.array(region[keys[k1]].historyOutputs['U1'].data)
            u2 = np.array(region[keys[k1]].historyOutputs['U2'].data)
        else:
            u1 = np.array([All_Tts[0, :, kk], All_Tts[1, :, kk]]).T
            u2 = np.array([All_Tts[0, :, kk], All_Tts[2, :, kk]]).T
        u = np.array([u1[:, 1], u2[:, 1]])
		
        # Get the node location
        node = np.where(measNodeCoords[:, 0] == nodeIdx)[0][0]
        
        # Angle
        phi = np.arcsin(r[0, node] / rmag[node])
        
        if abs(phi) > .41 and abs(phi) < .45:
            if rmag[node] > 40e-3:
                a = 1
            elif rmag[node] < 10e-3:
                a = 1
        
        # Spectrum
        time_L, time_S, fdr, fdp = getSpectrum(u, r[:,kk], rmag[kk], phi, time0, time1, spec_window, input_duration=input_duration, interp_factor=interp_factor)
        
        # Measured wave velocity
        try:
            meas_vL = (rmag[node] - refR) / (time_L - ref_tL)
            meas_vS = (rmag[node] - refR) / (time_S - ref_tS)
        except (FloatingPointError, ZeroDivisionError):
            meas_vL, meas_vS = 0, 0
        if np.isnan(meas_vL) or np.isnan(meas_vS):
            meas_vL, meas_vS = 0, 0
        
        # Undo wave propagation. Set to zero any values where the input spectrum
        # is close to zero (this is kills everything outside of approximate range
        # 2 MHz < f < 8 MHz).
        # Set to zero any values where the input spectrum is close to zero.
        inv_exp_prop_L = np.zeros((1, int(time_pts/2)), dtype=np.complex128)
        inv_exp_prop_S = np.zeros((1, int(time_pts/2)), dtype=np.complex128)
        for ii in range(idx1, idx2):
            inv_exp_prop_L[0, ii] = np.exp(+2.0j * np.pi * freq[ii] * time_L)
            inv_exp_prop_S[0, ii] = np.exp(+2.0j * np.pi * freq[ii] * time_S)
        inv_exp_prop_L = np.reshape(inv_exp_prop_L[0, idx1:idx2], (1, idx2-idx1))
        inv_exp_prop_S = np.reshape(inv_exp_prop_S[0, idx1:idx2], (1, idx2-idx1))
        
        # # Phase difference between propagation term e(iwt) and product of I^{-1} O
        # working_spec = np.dot(inv_in_freq_spec, fdr[idx1:idx2]) * np.sqrt(np.reshape(freq[idx1:idx2], (idx2-idx1, 1)))
        # working_phase = np.squeeze(np.angle(working_spec))
        # exp_phase = np.squeeze(np.angle(inv_exp_prop_L[0, :]))
    
        # phasediff = np.mod(exp_phase + working_phase + np.pi, 2*np.pi) - np.pi
        # phasediffstr = np.array_str(phasediff, precision=8, max_line_width=20*(idx2-idx1))
        # meanphasediff = nanmean(phasediff)
        
        # We have that (out_spec = in_spec * e(-iwt) * amps) and (amps = B * D * T)
        # T = 1 as we are in contact and are only looking at direct paths.
        # => D = e(iwt) * in_spec^{-1} * out_spec / B
        # where q = 3/4 for L; q = 5/4 for S.
        try:
            dir_L[kk] = np.dot(inv_exp_prop_L, np.dot(inv_in_freq_spec, fdr[idx1:idx2]) * np.sqrt(rmag[node]))[0][0]
            dir_S[kk] = np.dot(inv_exp_prop_S, np.dot(inv_in_freq_spec, fdp[idx1:idx2]) * np.sqrt(rmag[node]))[0][0]
        except (FloatingPointError, ZeroDivisionError):
            dir_L[kk], dir_S[kk] = 0, 0
        ####################################################################################################
        
        try:
            dL, dS = line_dir(phi, meas_vL/centre_freq, meas_vS/centre_freq), line_dir(phi, meas_vL/centre_freq, meas_vS/centre_freq, isLong=False)
        except (FloatingPointError, ZeroDivisionError):
            dL, dS = np.nan, np.nan
        
        # Write output
        f.write('{}, {}, {}, {}, {}, {}, {}, {}\n'.format(rmag[node], phi, dir_L[kk], dir_S[kk], meas_vL, meas_vS, dir_L[kk]-dL, dir_S[kk]-dS))#, phasediffstr))
        
        # Update progress bar based on platform.
        if "linux" in sys.platform and (kk+1)%(int(len(measNodeIdxs)/10)+1) == 0:
            updateProgress(t1, kk+1, len(measNodeIdxs))
        elif "win" in sys.platform and (kk+1)%(int(len(measNodeIdxs)/100)+1) == 0:
            updateProgress(t1, kk+1, len(measNodeIdxs))
    print('')

t2 = time.time()
print("\nTotal loop time: {}".format(secsToString(t2 - t0)))