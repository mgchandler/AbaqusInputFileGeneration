from odbAccess import openOdb
import sys

odb_name = sys.argv[1]
filename = odb_name[:-3]+'dat'

odb = openOdb(odb_name)

region = odb.steps['Step-1'].historyRegions
keys = region.keys()

with open(filename, 'w') as f:
	for kk in range(len(keys)):
		data = region[keys[kk]].historyOutputs['U2'].data
		filename = region[keys[kk]].name + '.dat'
		
		f.write(region[keys[kk]].name+'\n')
	
		for line in range(len(data)):
			f.write('{} {}\n'.format(data[line][0], data[line][1]))