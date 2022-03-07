# Writes .dat files containing timetraces for the specified .odb in the current
# working directory. Note that this script must be called from within Abaqus,
# i.e. at the command prompt, use command "abaqus python write_single_timetraces.py
# <odb_filename>.odb" in the directory containing the files.

# A timetrace is defined as U2 History Output data.

from odbAccess import openOdb
import sys

assert len(sys.argv) == 2

odb_name = sys.argv[1]
filename = odb_name[:-3]+'dat'

odb = openOdb(odb_name)

region = odb.steps['Step-1'].historyRegions
keys = region.keys()

with open(filename, 'w') as f:
	for kk in range(len(keys)):
		
		data_u1 = region[keys[kk]].historyOutputs['U1'].data
		data_u2 = region[keys[kk]].historyOutputs['U1'].data
		filename = region[keys[kk]].name + '.dat'
				
		# Write timetraces
		node_idx = filename.split('.')[1]
		f.write('Node '+node_idx+'\n')
	
		for line in range(len(data)):
			f.write('{} {} {}\n'.format(data_u1[line][0], data_u1[line][1], data_u2[line][1]))
			
		f.write('\n')