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
def getSpectrum(u, r, rmag, phi, time0, time1, spec_window, input_duration=1e-6, v=6317.012224890781, interp_factor=8, idx1=0, idx2=-1):
    # Get expected time of arrival of longitudinal and shear components
    sample_time = rmag / v + input_duration/2
		
    # Time-domain window to isolate signal of interest
    window = np.ones((u.shape[1], 1))
    window[time0 < sample_time - 1.05*input_duration/2] = 0
    window[time0 > sample_time + 1.05*input_duration/2] = 0
        
    # Radial and azimuthal displacement
    du = np.reshape(np.dot(u.T, r)/rmag, (u.shape[1], 1)) * window
        
    # Spectra
    fdu = np.zeros((du.shape[0], 1), dtype=complex)
    fdu[:int(u.shape[1]/2)+1] = 2 * np.fft.rfft(du, axis=0)
    fdu *= spec_window
        
    # Fourier interpolate to get higher resolution time of propagation
    idu = np.fft.ifft(np.append(interp_factor * fdu, np.zeros((fdu.shape[0]*(interp_factor-1), fdu.shape[1]))), axis=0)
		
    # Get the time of the maximum in the vicinity of the expected time of arrival, to compute actual
    # wave speed.
    time_of_arrival = time1[np.where(np.abs(idu) == 
     							np.max(np.abs(idu)))[0][0]] - input_duration/2
        
    return time_of_arrival, fdu[idx1:idx2]

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
        directivity = np.exp(3.0j/4.0 * np.pi) * cosT * (mu2 - 2 * sinT2) / F_0(sinT2, mu2)
    else:
        directivity = np.exp(5.0j/4.0 * np.pi) * mu2**(5/4) * np.sin(2*theta) * np.sqrt(mu2 * sinT2 - 1, dtype=np.complex128) / F_0(mu2*sinT2, mu2)
        
    return directivity

def F_0(z2, mu2):
    return (2 * z2 - mu2)**2 - 4 * z2 * np.sqrt(z2 - 1, dtype=np.complex128) * np.sqrt(z2 - mu2, dtype=np.complex128)



#%% Start. Determine if we are running from command line or not.

try:
    assert len(sys.argv) == 2
    odb_name = sys.argv[1]
    filename = odb_name[:-3]+'dir_k.dat'
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
c11 = modulus * (1. - poisson) / ((1. + poisson) * (1. - 2. * poisson))
c44 = modulus / (2. * (1. + poisson))



#%% Read input file to get coordinate locations.

if ".odb" in odb_name:
    odb = openOdb(odb_name, readOnly=True)
        
    numMeasNodes = len(odb.rootAssembly.instances['PART-1-1'].nodeSets['MEASURESET'].nodes)
    measNodeIdxs = np.zeros((numMeasNodes), dtype=int)
    measNodeCoords1 = np.zeros((numMeasNodes, 3), dtype=float)
    for node in range(numMeasNodes):
        measNodeIdxs[node] = int(odb.rootAssembly.instances['PART-1-1'].nodeSets['MEASURESET'].nodes[node].label)
        measNodeCoords1[node] = odb.rootAssembly.instances['PART-1-1'].nodeSets['MEASURESET'].nodes[node].coordinates
        
    measNodeCoords = np.array([measNodeCoords1[:,0], measNodeCoords1[:,1]]).T
else:
    measNodeIdxs = []
    measNodeCoords = np.zeros((0, 2))
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
                    coords = np.reshape([float(line[1]), float(line[2])], (1,2))
                    measNodeCoords = np.append(measNodeCoords, coords, axis=0)

# Radial distances
r = np.transpose(measNodeCoords)
rmag = np.zeros(r.shape[1])
for ii in range(r.shape[1]):
    rmag[ii] = np.linalg.norm(r[:, ii])

minR = np.where(rmag == min(rmag[rmag!=0.0]))[0][0]



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
else:
    region = odb.steps['Step-1'].historyRegions
    keys = region.keys()
    
    # Get the input signal.
    if 'CF2' in region[keys[0]].historyOutputs.keys():
        amplitude = np.array(region[keys[0]].historyOutputs['CF2'].data)
        ignoreNode1 = True
        
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
freq = ((np.array(range(fft_pts))) * fstep)
sig_window = np.ones((time0.shape[0], 1))
sig_window[time0 > input_duration] = 0
spec_window = np.ones((freq.shape[0], 1))
spec_window[freq > 2*centre_freq] = 0

in_time_sig = np.reshape(amplitude[:, 1], (amplitude.shape[0], 1)) * sig_window
in_freq_spec = np.zeros((in_time_sig.shape[0], 1), dtype=complex)
in_freq_spec[:int(in_time_sig.shape[0]/2 + 1)] = 2 * np.fft.rfft(in_time_sig, axis=0)
in_freq_spec *= spec_window

# Find out where the spectrum is non-zero for computing directivities later.
nonZeros = np.where((np.abs(in_freq_spec) > max(np.abs(in_freq_spec))/100))[0]
idx1, idx2 = nonZeros[0], nonZeros[-1]
for ii in range(idx1, idx2):
    if ii not in nonZeros:
        print("Warning: non zero region of input spectrum between {:.3e}Hz and {:.3e}Hz is not contiguous".format(freq[idx1], freq[idx2]))

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
if ".odb" in odb_name:
    u1 = np.array(region[keys[minR]].historyOutputs['U1'].data)
    u2 = np.array(region[keys[minR]].historyOutputs['U2'].data)
else:
    raise NotImplementedError("Check that minR is valid here, based on whether node 1 is included earlier on.")
    u1 = np.array([All_Tts[0, :, minR], All_Tts[1, :, minR]]).T
    u2 = np.array([All_Tts[0, :, minR], All_Tts[2, :, minR]]).T

u = np.array([u1[:, 1], u2[:, 1]])
phi = np.arcsin(r[0, minR] / rmag[minR])
r_perp = np.dot(np.array([[0, -1], [1, 0]]), r[:,minR])
ref_tL, _ = getSpectrum(u, r[:,minR], rmag[minR], phi, time0, time1, spec_window, input_duration=input_duration, v=6317.0122248907810, interp_factor=interp_factor)
ref_tS, _ = getSpectrum(u, r_perp,    rmag[minR], phi, time0, time1, spec_window, input_duration=input_duration, v=3110.2818131859126, interp_factor=interp_factor)
refR = rmag[minR]



#%% Compute directivities

dir_L = np.zeros(measNodeIdxs.shape[0], dtype=np.complex128)
dir_S = np.zeros(measNodeIdxs.shape[0], dtype=np.complex128)

t1 = time.time()
print("Time to loop start: {}".format(secsToString(t1 - t0)))

with open(filename, 'w') as f:
    f.write('rmag, phi, C dL, C dS, vL, vS\n')
    # f.write('{}\n'.format(np.array_str(freq[idx1:idx2], precision=1, max_line_width=20*(idx2-idx1))))
    for kk in range(measNodeIdxs.shape[0]):

        # Get the node location. Use this to refer to anything involving MeasureSet, i.e. r and rmag
        nodeIdx = int(keys[kk].split('.')[1])
        if nodeIdx not in measNodeIdxs:
            continue
        node = np.where(measNodeIdxs == nodeIdx)[0][0]

        ################################# REMOVE INDENTATIONS WITHIN HERE ##################################
        if 'U1' not in region[keys[kk]].historyOutputs.keys() or 'U2' not in region[keys[kk]].historyOutputs.keys() or abs(rmag[node]) < 5e-16:
            continue

        # Displacements
        u = np.zeros((2, time_pts))
        if ".odb" in odb_name:
            u[0, :] = np.asarray(region[keys[kk]].historyOutputs['U1'].data)[:, 1]
            u[1, :] = np.asarray(region[keys[kk]].historyOutputs['U2'].data)[:, 1]
        else:
            u[0, :] = np.asarray([All_Tts[0, :, kk], All_Tts[1, :, kk]])
            u[1, :] = np.asarray([All_Tts[0, :, kk], All_Tts[2, :, kk]])

        # Angle
        phi = np.arcsin(r[0, node] / rmag[node])

        # Spectrum
        r_perp = np.dot(np.array([[0, -1], [1, 0]]), r[:,node])
        time_L, fdr = getSpectrum(u, r[:,node], rmag[node], phi, time0, time1, spec_window, input_duration=input_duration, v=6317.0122248907810, interp_factor=interp_factor, idx1=idx1, idx2=idx2)
        time_S, fdp = getSpectrum(u, r_perp,    rmag[node], phi, time0, time1, spec_window, input_duration=input_duration, v=3110.2818131859126, interp_factor=interp_factor, idx1=idx1, idx2=idx2)

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
        omega_t_L = 2.0 * np.pi * freq * time_L
        omega_t_S = 2.0 * np.pi * freq * time_S
        inv_exp_prop_L = np.zeros((1, idx2-idx1), dtype=np.complex128)
        inv_exp_prop_S = np.zeros((1, idx2-idx1), dtype=np.complex128)
        inv_exp_prop_L[0] = np.exp(omega_t_L[idx1:idx2] * 1.0j)
        inv_exp_prop_S[0] = np.exp(omega_t_S[idx1:idx2] * 1.0j)

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
        try:
            dir_L[node] = np.dot(inv_exp_prop_L, np.dot(inv_in_freq_spec, fdr) * np.sqrt(rmag[node]))[0][0] * c44 / np.sqrt(np.pi)# * np.sqrt(omega_t_L[idx1:idx2]))[0][0] * c11
            dir_S[node] = np.dot(inv_exp_prop_S, np.dot(inv_in_freq_spec, fdp) * np.sqrt(rmag[node]))[0][0] * c44 / np.sqrt(np.pi)# * np.sqrt(omega_t_S[idx1:idx2]))[0][0] * c44
            # dir_L[node] = np.dot(inv_exp_prop_L, np.dot(inv_in_freq_spec, fdr) * np.sqrt(omega_t_L[idx1:idx2]))[0][0] * c11
            # dir_S[node] = np.dot(inv_exp_prop_S, np.dot(inv_in_freq_spec, fdp) * np.sqrt(omega_t_S[idx1:idx2]))[0][0] * c44
        except (FloatingPointError, ZeroDivisionError):
            dir_L[node], dir_S[node] = 0, 0
        ####################################################################################################
        
        try:
            dL, dS = line_dir(phi, meas_vL/centre_freq, meas_vS/centre_freq), line_dir(phi, meas_vL/centre_freq, meas_vS/centre_freq, isLong=False)
        except (FloatingPointError, ZeroDivisionError):
            dL, dS = np.nan, np.nan
        
        # Write output
        f.write('{}, {}, {}, {}, {}, {}\n'.format(rmag[node], phi, dir_L[node], dir_S[node], meas_vL, meas_vS))#, phasediffstr))
        
        # Update progress bar based on platform.
        if "linux" in sys.platform and (kk+1)%(int(len(measNodeIdxs)/10)+1) == 0:
            updateProgress(t1, kk+1, len(measNodeIdxs))
        elif "win" in sys.platform and (kk+1)%(int(len(measNodeIdxs)/100)+1) == 0:
            updateProgress(t1, kk+1, len(measNodeIdxs))
        
    print('')

t2 = time.time()
print("\nTotal loop time: {}".format(secsToString(t2 - t0)))