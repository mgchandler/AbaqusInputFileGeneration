# Writes .dat files containing timetraces for the specified .odb in the current
# working directory. Note that this script must be called from within Abaqus,
# i.e. at the command prompt, use command "abaqus python write_single_timetraces.py
# <odb_filename>.odb <output_field>" in the directory containing the files.

# A timetrace is defined as each node's History Output data. Valid output fields are
# U1 (x-displacement), U2 (z-displacement), CF2 (z-force, concentrated). Timetraces
# will be written for all nodes where there is data for the specified output field.

from __future__ import print_function # For progress bar
import numpy as np
from odbAccess import openOdb
from operator import attrgetter, itemgetter
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



allowed_fields = ['U1', 'U2', 'CF2']

assert len(sys.argv) > 1
if len(sys.argv) >= 3:
    output_fields = sys.argv[2:]
    for output_field in output_fields:
        if output_field not in allowed_fields:
            raise ValueError("Invalid output field {}: choose from {}".format(output_field, ", ".join(allowed_fields)))
else:
	output_fields = ['U2']

odb_name = sys.argv[1]
filename = odb_name[:-4]+'_'+'_'.join(output_fields)+'.dat'

odb = openOdb(odb_name, readOnly=True)

region = odb.steps['Step-1'].historyRegions
keys = region.keys()

sets = odb.rootAssembly.instances['PART-1-1'].nodeSets
num_els = sum(np.char.count(sets.keys(), 'EL'))
nodes_per_el = len(sets['EL1'].nodes)
nodes_in_el = np.zeros((num_els, nodes_per_el), dtype=int)
for el in range(num_els):
	for node in range(nodes_per_el):
		nodes_in_el[el, node] = int(sets['EL{}'.format(el+1)].nodes[node].label)

doMeasureSet = False
if "MEASURESET" in sets.keys():
	doMeasureSet = True
	num_nodes = len(sets["MEASURESET"].nodes)
	nodes_in_measureset = np.zeros((num_nodes), dtype=int)
	for node in range(num_nodes):
		nodes_in_measureset[node] = int(sets["MEASURESET"].nodes[node].label)

t1 = time.time()
with open(filename, 'w') as f:
    for kk in range(len(keys)):
        # Get only history outputs of interest. Skip this node if they aren't all available.
        try:
            data = itemgetter(*output_fields)(region[keys[kk]].historyOutputs)
        except KeyError:
            continue
        
        # When only one output requested, itemgetter returns the item only, and a tuple otherwise.
        # Ensure that we always have a tuple.
        if not isinstance(data, tuple):
            data = (data,)
        
        # Get the history data only.
        data = list(map(attrgetter("data"), data))
        
        # Write timetraces
        node_idx = region[keys[kk]].name.split('.')[1]
        # If we have node sets of transducer elements.
        if int(node_idx) in nodes_in_el:
            el_idx = np.where(nodes_in_el == int(node_idx))[0][0] + 1
            f.write('Node-'+node_idx+' El{}\n'.format(el_idx))
            
            for line in range(len(data[0])):
                field_vals = []
                for output_field in range(len(data)):
                    field_vals.append(data[output_field][line][1])
                f.write('{} {}\n'.format(data[0][line][0], ' '.join(field_vals)))
                
            f.write('\n')
        # If we have a MeasureSet.
        elif doMeasureSet:
            if int(node_idx) in nodes_in_measureset:
                set_idx = np.where(nodes_in_measureset == int(node_idx))[0][0]
                f.write('{}.Node-'.format(kk+1)+node_idx+' MeasureSet coords=({:.8f},{:.8f},{:.8f})\n'.format(sets["MEASURESET"].nodes[set_idx].coordinates[0], sets["MEASURESET"].nodes[set_idx].coordinates[1], sets["MEASURESET"].nodes[set_idx].coordinates[2]))
                
                for line in range(len(data[0])):
                    field_vals = []
                    for output_field in range(len(data)):
                        field_vals.append('{}'.format(data[output_field][line][1]))
                    f.write('{} {}\n'.format(data[0][line][0], ' '.join(field_vals)))
                    
                f.write('\n')
        
        if (kk+1)%int(len(keys)/10) == 0:
            updateProgress(t1, kk, len(keys))