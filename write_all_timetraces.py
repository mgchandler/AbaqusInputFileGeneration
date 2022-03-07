# Writes .dat files containing timetraces for all .odb files in the current
# working directory. Note that this script must be called from within Abaqus,
# i.e. at the command prompt, use command "abaqus python write_all_timetraces.py"
# in the directory containing the files.

# A timetrace is defined as U2 History Output data.

from odbAccess import openOdb
import os

path = os.getcwd()
for odb_name in os.listdir('.'):
# Get all odb files in the current working directory
	if odb_name.endswith('odb'):
		data_filename = odb_name[:-3] + 'dat'
		odb = openOdb(odb_name)
		
		# Get complete set of history data from this transmission.
		region = odb.steps['Step-1'].historyRegions
		keys = region.keys()
		
		# Start writing
		with open(data_filename, 'w') as f:
			for kk in range(len(keys)):
				# Get each receiver
				data = region[keys[kk]].historyOutputs['U2'].data
				filename = region[keys[kk]].name + '.dat'
				
				# Write timetraces
				node_idx = filename.split('.')[1]
				f.write('Node '+node_idx+'\n')
			
				for line in range(len(data)):
					f.write('{} {}\n'.format(data[line][0], data[line][1]))
				f.write('\n')
		print("Data from file '{}.odb' written".format(data_filename))