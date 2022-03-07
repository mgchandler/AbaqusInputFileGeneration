# coding=utf_8

# Writes .dat files containing directivities for the specified .odb in the current
# working directory. Note that this script must be called from within Abaqus,
# i.e. at the command prompt, use command "abaqus python write_all_direcitivities.py
# <odb_filename>.odb" in the directory containing the files.

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
    
    
    
#%% Start. Determine if we are running from command line or not.

try:
    assert len(sys.argv) == 2
    odb_name = sys.argv[1]
    isCommand = True
    from odbAccess import openOdb
except AssertionError:
    odb_name = 'D_rdm1_30npw.dat'
    os.chdir(r"C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\Abaqus\FMC Generation\v9\Output\Phase Difference Study\Directivity and Speed convergence")
    isCommand = False
    import matplotlib.pyplot as plt

t0 = time.time()

filename = odb_name[:-3]+'phasediff.dat'
inpname = odb_name[:-3]+'inp'
datname = odb_name[:-3]+'tt.dat'

# Input parameters
v_L = 6292.054788836226 #6263.0
v_S = 3063.23298960871 #3083.0
input_duration = 1e-6
centre_freq = 5e6
no_cycles = 5.0
interp_factor = 4



#%% Read input file to get coordinate locations, and read in time traces if running in IDE.

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
perp_matrix = np.array([[0, -1], [1, 0]])
r_perp = np.dot(perp_matrix, r)
rmag = np.zeros(r.shape[1])
for ii in range(r.shape[1]):
    rmag[ii] = np.linalg.norm(r[:, ii])



#%% Read timetraces and get directivities

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
            elif count%2 == 1:
                if abs(rmag[kk] - float(line[0])) < 1e-8:
                    All_Tts[0, :, kk] = timevec
                    All_Tts[1, :, kk] = np.array([float(val) for val in line[1:]])
                else:
                    continue
            else:
                if (kk != 0 and abs(np.arcsin(r[0, kk] / rmag[kk]) - float(line[0])) < 1e-8) or abs(float(line[0])) < 2e-13:
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

else:
    odb = openOdb(odb_name, readOnly=True)
    
    region = odb.steps['Step-1'].historyRegions
    keys = region.keys()

if 1 in measNodeIdxs:
    readInput = True
    firstNode = False
else:
    readInput = False
    firstNode = True
	
t1 = time.time()
print("Time to loop start: {}".format(secsToString(t1 - t0)))

with open(filename, 'w') as f:
    f.write('rmag, phi, C dL, C dS, vL, vS, DPhase(f) : f = ')
    for kk in range(measNodeIdxs.shape[0]):
        
        ################################# REMOVE INDENTATIONS WITHIN HERE ##################################
        nodeIdx = measNodeIdxs[kk]
        
        # if rmag[kk] > 6e-3:
        #     continue
        
        # Displacements
        if ".odb" in odb_name:
            try:
                assert int(keys[kk].split('.')[1]) == nodeIdx 
            except AssertionError:
                print("Input file node {} does not match odb file node {}".format(keys[kk].split('.')[1], nodeIdx))
                raise
            u1 = np.array(region[keys[kk]].historyOutputs['U1'].data)
            u2 = np.array(region[keys[kk]].historyOutputs['U2'].data)
        else:
            u1 = np.array([All_Tts[0, :, kk], All_Tts[1, :, kk]]).T
            u2 = np.array([All_Tts[0, :, kk], All_Tts[2, :, kk]]).T
        
        u = np.array([u1[:, 1], u2[:, 1]])
        umag = np.zeros(u.shape[1])
        for ii in range(u.shape[1]):
            umag[ii] = np.linalg.norm(u[:, ii])
        
        # If we're reading the input, the first node will contain the input signal. Read it to obtain fdz
        # and then continue to the next node.
        if readInput:
            print('Using first node')
            assert(abs(rmag[kk]) < 2e-13)
            time0 = np.array(u1[:, 0])
            time_pts = u.shape[1]
            time_step = time0[1] - time0[0]
            
            fft_pts = time_pts
            fstep = 1.0/(time0[-1])
            freq = ((np.array(range(int(fft_pts/2)))) * fstep)
            idx1 = np.where(np.abs(freq - 2e6) == np.min(np.abs(freq - 2e6)))[0][0]-1
            idx2 = np.where(np.abs(freq - 8e6) == np.min(np.abs(freq - 8e6)))[0][0]+1
            spec_window = np.ones((freq.shape[0], 1))
            spec_window[freq > 2*centre_freq] = 0

            in_time_sig = np.reshape(u[1, :], (umag.shape[0], 1))
            in_freq_spec = np.fft.fft(in_time_sig, axis=0)
            in_freq_spec = in_freq_spec[:int(in_time_sig.shape[0] / 2)] * spec_window
            inv_in_freq_spec = np.diag(1.0/np.squeeze(in_freq_spec[idx1+1:idx2+1]))
            
            interp_freq_spec = np.zeros((in_freq_spec.shape[0]*interp_factor, 1), dtype=np.complex128)
            interp_freq_spec[:in_freq_spec.shape[0]] = interp_factor * in_freq_spec
            # interp_freq_spec[-in_freq_spec.shape[0]:] = np.flip(in_freq_spec)
            idu = np.fft.ifft(np.append(interp_factor * in_freq_spec, np.zeros((in_freq_spec.shape[0]*(interp_factor-1), in_freq_spec.shape[1]))), axis=0)
            input_duration = 2 * np.linspace(time0[0], time0[-1], idu.shape[0])[np.where(np.abs(idu) == np.max(np.abs(idu)))[0][0]]
            
            time1 = np.linspace(time0[0], time0[-1], idu.shape[0])#time0[np.array(range(len(time0)))%2 == 0]
            # time1 = time1[:len(in_freq_spec)]
            f.write('{}\n'.format(np.array_str(freq[idx1:idx2], precision=1, max_line_width=20*(idx2-idx1))))
            
            readInput = False
            
            continue
            
        # Create the input signal. Do this only once, if we are on the first node.
        elif firstNode:
            firstNode = False
            
            print("No first node")
            
            # Make sure the input signal has the same shape as the measured signals
            # we will be working with later on.
            time_pts = u.shape[1]
            time_step = u1[-1, 0] / (time_pts-1)
            time = np.array(range(time_pts)) * time_step
            tmax = max(time)
            tmid = tmax/2
            half_width_fract = no_cycles / (centre_freq * tmax * 2.0)
            input_duration = no_cycles / centre_freq
            
            carrier = np.sin(2.0*np.pi * centre_freq * (time - tmid))
            window = np.exp(-((np.linspace(0, 1, time_pts) - .5) / (half_width_fract / np.sqrt(-np.log(10**(-40/20))))) ** 2.0)
            time_sig = carrier * window
            
            fft_pts = time_pts
            fstep = 1.0/(fft_pts * time_step)
            freq = ((np.array(range(fft_pts/2))) * fstep)
            in_freq_spec = np.fft.rfft(time_sig, axis=0)
            in_freq_spec = in_freq_spec[:int(fft_pts / 2)]
            in_freq_spec = in_freq_spec * np.exp(2j * np.pi * freq * tmid)
            in_time_sig = np.real(np.fft.ifft(in_freq_spec, fft_pts)) * 2.0
            in_time_sig = in_time_sig[:time_pts]
            sf = 1.0 / np.max(np.abs(in_time_sig))
            in_time_sig = in_time_sig * sf
            in_freq_spec = in_freq_spec * sf
            
            idx1 = np.where(np.abs(freq - 2e6) == np.min(np.abs(freq - 2e6)))[0][0]
            idx2 = np.where(np.abs(freq - 8e6) == np.min(np.abs(freq - 2e6)))[0][0]
            f.write('{}\n'.format(freq[idx1:idx2]))
            
            inv_in_freq_spec = np.diag(1/in_freq_spec[idx1:idx2])
		
        # Get the node location
        node = np.where(measNodeCoords[:, 0] == nodeIdx)[0][0]
        
        # Angles
        phi = np.arcsin(r[0, node] / rmag[node])
        
        # # Get trig denominator. Do this separately as there are cases where b == 0.
        # b = np.reshape((umag * rmag[node]), (umag.shape[0], 1))
        # # We'll define dr, dp := 0 where umag == 0. This is obtained in the numerator, 
        # # so set denominator to non-zero to avoid NaNs.
        # b[b == 0] = 1.0
        # cosarg = -np.dot(u.T, np.reshape(r[:, node], (2,1))) / b
        # sinarg = np.dot(u.T, np.reshape(r_perp[:, node], (2,1))) / b
    
        # Get expected time of arrival of longitudinal and shear components
        long_sample_time = rmag[kk] / v_L + input_duration/2
        shear_sample_time = rmag[kk] / v_S + input_duration/2
		
        l_window = np.ones((umag.shape[0], 1))
        l_window[time0 < long_sample_time - 1.05*input_duration/2] = 0
        l_window[time0 > long_sample_time + 1.05*input_duration/2] = 0
        s_window = np.ones((umag.shape[0], 1))
        s_window[time0 < shear_sample_time - 1.05*input_duration/2] = 0
        s_window[time0 > shear_sample_time + 1.05*input_duration/2] = 0
        
        dr = np.reshape(np.dot(-u.T, r[:,kk])/rmag[kk], (u.shape[1], 1)) * l_window
        dp = np.reshape(np.dot(u.T, r_perp[:,kk])/rmag[kk], (u.shape[1], 1)) * s_window
        
        fdr = np.fft.rfft(dr, axis=0)
        fdp = np.fft.rfft(dp, axis=0)
        # fdz = np.fft.rfft(dz, axis=0)
        fdr = fdr[:int(dr.shape[0] / 2)] * spec_window
        fdp = fdp[:int(dp.shape[0] / 2)] * spec_window
        # fdz = fdz[:int(dz.shape[0] / 2)]
        
        idr = np.fft.ifft(np.append(interp_factor * fdr, np.zeros((fdr.shape[0]*(interp_factor-1), fdr.shape[1]))), axis=0)
        idp = np.fft.ifft(np.append(interp_factor * fdp, np.zeros((fdp.shape[0]*(interp_factor-1), fdp.shape[1]))), axis=0)
        # dz = np.fft.ifft(fdz, axis=0)
		
        # Get the time of the maximum in the vicinity of the expected time of arrival, to compute actual
        # wave speed.
        time_L = time1[np.where(np.abs(idr) == 
         							np.max(np.abs(idr)))[0][0]] - input_duration/2
        time_S = time1[np.where(np.abs(idp) == 
        							 np.max(np.abs(idp)))[0][0]] - input_duration/2
        meas_vL = rmag[node] / time_L
        meas_vS = rmag[node] / time_S
        
        # Set to zero any values where the input spectrum is close to zero. This really needs to be more
        # robust, but by inspection idxs 34:133 catch all frequencies 2MHz < f < 8MHz for this data.
        inv_exp_prop_L = np.zeros((1, int(time_pts/2)), dtype=np.complex128)
        inv_exp_prop_S = np.zeros((1, int(time_pts/2)), dtype=np.complex128)
        for ii in range(idx1, idx2):
            inv_exp_prop_L[0, ii] = np.exp(2.0j * np.pi * freq[ii] * time_L)
            inv_exp_prop_S[0, ii] = np.exp(2.0j * np.pi * freq[ii] * time_S)
        inv_exp_prop_L = np.reshape(inv_exp_prop_L[0, idx1:idx2], (1, idx2-idx1))
        inv_exp_prop_S = np.reshape(inv_exp_prop_S[0, idx1:idx2], (1, idx2-idx1))
    
        working_spec = np.dot(inv_in_freq_spec, np.reshape(fdr[idx1:idx2], (fdr[idx1:idx2].shape[0], 1)))# * np.sqrt(np.reshape(freq[1:], (fdr.shape[0]-1, 1)))
        # working_freq = freq[idx1:idx2]
        working_phase = np.squeeze(np.angle(working_spec))
        exp_phase = np.squeeze(np.angle(inv_exp_prop_L[0, :]))
    
        phasediff = np.mod(exp_phase + working_phase + np.pi, 2*np.pi) - np.pi
        phasediffstr = np.array_str(phasediff, precision=8, max_line_width=20*(idx2-idx1))
        meanphasediff = nanmean(phasediff)
        
        # We have that (out_spec = in_spec * e(-iwt) * amps) and (amps = B * D * T)
        # T = 1 as we are in contact and are only looking at direct paths.
        # => D = e(iwt) * in_spec^{-1} * out_spec / B
        dir_L_phased = np.dot(np.dot(inv_exp_prop_L, inv_in_freq_spec), fdr[idx1:idx2])[0][0] * np.sqrt(rmag[node])
        dir_S_phased = np.dot(np.dot(inv_exp_prop_S, inv_in_freq_spec), fdp[idx1:idx2])[0][0] * np.sqrt(rmag[node])
		####################################################################################################
        
        # Write output
        f.write('{}, {}, {}, {}, {}, {}, {}\n'.format(rmag[node], phi, dir_L_phased, dir_S_phased, meas_vL, meas_vS, phasediffstr))
        
        # Update progress bar based on platform.
        if "linux" in sys.platform and (kk+1)%(int(len(measNodeIdxs)/10)+1) == 0:
            updateProgress(t1, kk+1, len(measNodeIdxs))
        elif "win" in sys.platform and (kk+1)%(int(len(measNodeIdxs)/100)+1) == 0:
            updateProgress(t1, kk+1, len(measNodeIdxs))

t2 = time.time()
print("\nTotal runtime: {}".format(secsToString(t2 - t0)))