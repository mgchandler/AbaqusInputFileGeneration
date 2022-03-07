# Writes a .dat file containing the field displacement magnitude for the specified 
# .odb in the current working directory, at a specified frame. Note that this script
# must be called from within Abaqus, i.e. at the command prompt, use command "abaqus
# python write_single_frame.py <odb_filename>.odb <frame_no>" in the directory
# containing the files.

from odbAccess import openOdb
import sys

assert len(sys.argv) == 3

odb_name = sys.argv[1]
frame_no = int(sys.argv[2])
filename = odb_name[:-3]+'dat'

odb = openOdb(odb_name)

frame = odb.steps['Step-1'].frames[frame_no]

with open(filename, 'w') as f:
	
	data = odb.steps['Step-1'].frames[frame_no].fieldOutputs['U'].values
	all_nodes = odb.rootAssembly.instances['Part-1-1'].nodes
	
	f.write('NodeID    x    z    U1    U2    UMag\n')
	for node in range(len(data)):
		this_U = data[node]
		this_node = all_nodes[node]
		f.write('{}    {}    {}    {}    {}    {}\n'.format(
			this_U.nodeLabel, this_node.coordinates[0], this_node.coordinates[2], 
			this_U.data[0], this_U.data[1], this_U.magnitude))