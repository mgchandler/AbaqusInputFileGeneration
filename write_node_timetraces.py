# Writes .dat files containing timetraces for the specified .odb in the current
# working directory. Note that this script must be called from within Abaqus,
# i.e. at the command prompt, use command "abaqus python write_single_timetraces.py
# <odb_filename>.odb" in the directory containing the files.

# A timetrace is defined as U2 History Output data.

from __future__ import print_function # For progress bar
import numpy as np
from odbAccess import openOdb
import sys
import time

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

def vecToString(vector):
	vector = np.squeeze(vector)
	assert len(vector.shape) == 1
	assert vector.shape[0] > 1
	string = "{:.8g}".format(vector[0])
	for value in vector[1:]:
		string += ", {:.8g}".format(value)
	return string

assert len(sys.argv) == 2
input_duration = 1e-6
centre_freq = 5e6
interp_factor = 4

odb_name = sys.argv[1]
inpname = odb_name[:-3]+'inp'
filename = odb_name[:-4]+'.tt.dat'

odb = openOdb(odb_name, readOnly=True)

region = odb.steps['Step-1'].historyRegions
keys = region.keys()

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
# measNodeIdxs = measNodeIdxs[:6]

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
    
    time0 = np.array(region[keys[0]].historyOutputs['U1'].data)[:, 0]
    amplitude = np.zeros((time0.shape[0], 2))
    amplitude[:,0] = time0
    amplitude[:,1] = np.interp(time0, readAmplitude[:,0], readAmplitude[:,1], right=0)

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

nonZeros = np.where((np.abs(in_freq_spec) > max(np.abs(in_freq_spec))/100))[0]
idx1, idx2 = nonZeros[0], nonZeros[-1]

inv_in_freq_spec = np.diag(1.0/np.squeeze(in_freq_spec[idx1+1:idx2+1]))

interp_freq_spec = np.zeros((in_freq_spec.shape[0]*interp_factor, 1), dtype=np.complex128)
interp_freq_spec[:in_freq_spec.shape[0]] = interp_factor * in_freq_spec
idu = np.fft.ifft(np.append(interp_factor * in_freq_spec, np.zeros((in_freq_spec.shape[0]*(interp_factor-1), in_freq_spec.shape[1]))), axis=0)
input_duration = 2 * np.linspace(time0[0], time0[-1], idu.shape[0])[np.where(np.abs(idu) == np.max(np.abs(idu)))[0][0]]
        
time1 = np.linspace(time0[0], time0[-1], idu.shape[0])

t1 = time.time()
print("Starting file write")
with open(filename, 'w') as f:
    outName = region[keys[0]].historyOutputs.keys()[0]
    if 'CF2' in outName:
        ignoreNode0 = True
    else:
        ignoreNode0 = False
    timestr = vecToString(time0)
    f.write('t, {}\n'.format(timestr))
    f.write('i, {}\n'.format(vecToString(in_time_sig[:,0])))
    for kk in range(len(measNodeIdxs)):
        if ignoreNode0:
            k1 = kk + 1
        else:
            k1 = kk
        node_idx = region[keys[k1]].name.split('.')[1]
        assert int(node_idx) == int(measNodeCoords[kk, 0])
        		
        data_u1 = np.array(region[keys[k1]].historyOutputs['U1'].data)
        data_u1 = vecToString(data_u1[:, 1])
        data_u2 = np.array(region[keys[k1]].historyOutputs['U2'].data)
        data_u2 = vecToString(data_u2[:, 1])
          		
        # Write timetraces
        try:
            phi = np.arcsin(r[0, kk] / rmag[kk])
        except FloatingPointError:
            phi = 0.0
          		
        f.write("{:.8f}, {}\n".format(rmag[kk], data_u1))
        f.write("{:.8f}, {}\n".format(phi, data_u2))
        
        # Update progress bar based on platform.
        if "linux" in sys.platform and (kk+1)%(int(len(measNodeIdxs)/10)+1) == 0:
            updateProgress(t1, kk, len(measNodeIdxs))
        elif "win" in sys.platform and (kk+1)%(int(len(measNodeIdxs)/1000)+1) == 0:
            updateProgress(t1, kk, len(measNodeIdxs))
    print('')