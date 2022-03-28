# Integrates across all history outputs contained within a single file to produce
# FMC data in separate files. Call with command "python compile_timetraces.py
# <filename>" or "python compile_timetraces.py <filename>.dat". <filename> is
# assumed to have form "<jobname>_16" or "<jobname>_16_U2" to specify that this is
# transmitting element 16, containing displacement in the y-axis. Files are expected
# to have format 
# --- <filename.dat> ---
# L1    : Node-<node1> El<el1>
# L2    : time0 value0
# L3    : time1 value1
# ...     ...
# LI+2  : timeI valueI
# LI+3  :
# LI+4  : Node-<node2> El<el2>
# LI+5  : time0 value0
# ...     ...
# L2I+5 : timeI valueI
# ...     ...
# --- <end of file> ---
# where <node> is the node number from the .odb or .inp file, and El<el> is the set
# which it is contained by.
#
# It is assumed that all nodes have the same time axis, although this assumption is
# not checked. It is also assumed that nodes within each element are equally spaced,
# however this assumption is also not checked.

import numpy as np
import sys

Jobname = sys.argv[1]
try:
    tx = int(Jobname.split('_')[-1])
except ValueError:
    tx = int(Jobname.split('_')[-2])

if Jobname[-4:] != ".dat":
    filename = Jobname + ".dat"
else:
    filename = Jobname

# Loop through file once to get timetrace length, number of nodes, elements associated with each node and the time axis.
el_labels = []
num_nodes = 0
num_time_pts = 0
time_counter = 0
with open(filename, 'r') as f:
    for count, line in enumerate(f):
        # Start of a new node.
        if "El" in line:
            num_nodes += 1
            el_labels.append(int(line.split(' ')[1][2:]))
            if num_nodes == 2:
                time_axis = np.zeros((num_time_pts, 1))
        # First node: work out length of timetrace.
        elif num_nodes == 1 and line != '\n':
            num_time_pts += 1
        # Second node: now we know tt length, get the time axis.
        elif num_nodes == 2 and line != '\n':
            time_axis[time_counter, 0] = np.double(line.split(' ')[0])
            time_counter += 1

num_els = np.unique(el_labels).shape[0]

# Loop through a second time to get the actual data. One loop would require appending; two enables pre-defining array size.
# For larger files, this is much more efficient; for smaller files, appending would probably be quicker but the actual runtime difference would be negligible.
timetraces = np.zeros((2, num_time_pts, num_els))
timetraces[0, :, :] = time_axis
node_counter = -1
with open(filename, 'r') as f:
    for line in f:
        # New node: reset for next column in array.
        if "Node" in line:
            node_counter += 1
            time_counter = 0
            continue
        # Integrate all nodes across the length of the array. Note that nodes within the array are assumed to have equal spacing.
        elif line != '\n':
            timetraces[1, time_counter, el_labels[node_counter]-1] += np.double(line.split(' ')[1])
            time_counter += 1
        else:
            continue

# Remove any trend in the data.
for el in range(num_els):
    poly = np.polynomial.Polynomial.fit(timetraces[0, :, el], timetraces[1, :, el], 1).convert().coef
    timetraces[1, :, el] -= poly[1] * timetraces[0, :, el] + poly[0]

# Write out the FMC data.
for rx in range(num_els):
    with open('tx{}-rx{}.dat'.format(tx, rx+1), 'w') as file:
        for line in range(timetraces.shape[1]):
            file.write('{} {}\n'.format(timetraces[0, line, rx], timetraces[1, line, rx]))