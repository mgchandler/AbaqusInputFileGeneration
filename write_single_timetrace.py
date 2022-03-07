# Writes .dat files containing timetraces for the specified .odb in the current
# working directory. Note that this script must be called from within Abaqus,
# i.e. at the command prompt, use command "abaqus python write_single_timetraces.py
# <odb_filename>.odb <output_field>" in the directory containing the files.

# A timetrace is defined as each node's History Output data. Valid output fields are
# U1 (x-displacement), U2 (z-displacement), CF2 (z-force, concentrated). Timetraces
# will be written for all nodes where there is data for the specified output field.

from odbAccess import openOdb
import sys

output_fields = ['U1', 'U2', 'CF2']

assert len(sys.argv) > 1
if len(sys.argv) == 3:
	output_field = sys.argv[2]
	if output_field not in output_fields:
		raise ValueError("Invalid output field {}: choose from {}".format(output_field, ", ".join(output_fields)))
else:
	output_field = 'U2'

odb_name = sys.argv[1]
filename = odb_name[:-4]+'_'+output_field+'.dat'

odb = openOdb(odb_name, readOnly=True)

region = odb.steps['Step-1'].historyRegions
keys = region.keys()

with open(filename, 'w') as f:
	for kk in range(len(keys)):
		if output_field in region[keys[kk]].historyOutputs.keys():
			data = region[keys[kk]].historyOutputs[output_field].data
			filename = region[keys[kk]].name + '.dat'
			
			# Write timetraces
			node_idx = filename.split('.')[1]
			f.write('Node '+node_idx+'\n')
			
			for line in range(len(data)):
				f.write('{} {}\n'.format(data[line][0], data[line][1]))
			
			f.write('\n')