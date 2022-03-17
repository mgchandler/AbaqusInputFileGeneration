# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 16:02:09 2022

@author: mc16535
    
Call this script from the command line with "python plot_timetraces.py <filename> 
**kwargs". Plots data contained within the file <filename>, or any subset of this
described by **kwargs.

Parameters
----------
filename : string
    Name of the file containing data. Assumed to be in the working directory.
    
node_start : int
    (default 0) The first node in the range of nodes contained in the file which will be
    plotted. This specifies the node as it appears in the file, not the node
    label (i.e. the file may record data for nodes 3000-3010: node 0 refers to
    3000, 1 to 3001, etc.) Note that only consecutive nodes can be plotted in
    one figure.
    
node_end : int
    (default -1) The last node in the range to be plotted.
    
col_start : int
    (default 1) The first column in the range to be plotted. Note that only consecutive
    ranges can be plotted.
    
col_end : int
    (default -1) The last column to be plotted.

Notes
----------
Files are assumed to be formatted in the following ways:

--- <filename.dat> ---
L2    : time0, avalue0, bvalue0, ...
L3    : time1, avalue1, bvalue1, ...
...     ...
LI    : timeI, avalueI, bvalueI, ...
--- <end of file> ---

OR

--- <filename.dat> ---
L1    : Node 1
L2    : time0, avalue0, bvalue0, ...
L3    : time1, avalue1, bvalue1, ...
...     ...
LI+2  : timeI, avalueI, bvalueI, ...
LI+3  :
LI+4  : Node 2
LI+5  : time0, avalue0, bvalue0, ...
...     ...
L2I+5 : timeI, avalueI, bvalueI, ...
...     ...
--- <end of file> ---

OR

--- <filename.dat> ---
L1    : Node 1
L2    : column1, column2, column3, ...
L3    : time0, avalue0, bvalue0, ...
L4    : time1, avalue1, bvalue1, ...
...     ...
LI+3  : timeI, avalueI, bvalueI, ...
LI+4  :
LI+5  : Node 2
LI+6  : column1, column2, column3, ...
LI+7  : time0, avalue0, bvalue0, ...
...     ...
L2I+7 : timeI, avalueI, bvalueI, ...
...     ...
--- <end of file> ---

"""

import matplotlib.pyplot as plt
import numpy as np
import sys

DELIMITERS = [',', ';', '|', '\t', '  ', ' ']

#%%
def main(filename, node_start=0, node_end=-1, col_start=1, col_end=-1):
    node_start = int(node_start)
    node_end = int(node_end)
    col_start = int(col_start)
    col_end = int(col_end)
    
    if node_end != -1 and node_end < node_start:
        raise ValueError("Last node index comes before first node index.")
    if col_end != -1 and col_end < col_start:
        raise ValueError("Last column index comes before first column index.")
    if node_start < 0 or col_start < 0:
        raise ValueError("Start index too small.")
    
    
    
    #%% Work out what the delimiter is. Expect line 3 to always contain data.
    with open(filename, 'r') as file:
        line3 = file.readlines()[2]
        delimiter = None
        for d in DELIMITERS:
            if d in line3:
                delimiter = d
                if delimiter == ' ':
                    print("Warning - delimiter==' ', labels may break.")
                
                break
        
        num_cols = len(line3.split(delimiter))
    
    
    
    #%% Count the number of nodes.
    num_nodes = 1
    with open(filename, 'r') as file:
        for count, line in enumerate(file):
            if "Node" in line and count != 0:
                num_nodes += 1
    
    if node_end == -1:
        node_end = num_nodes
    if col_end == -1:
        col_end = num_cols
    
    if node_end > num_nodes or col_end > num_cols:
        raise ValueError("End index too large.")
    
    
    
    #%% Read in the data.
    data = np.full((num_nodes, int(count/num_nodes)+1, num_cols), np.nan)
    labels = np.full((num_nodes, num_cols), "")
    node = -1
    ii = 0
    with open(filename, 'r') as file:
        label_prefix = ""
        for line in file:
            if "Node" in line:
                label_prefix = line
                node += 1
                ii = 0
                continue
            
            if line == "\n":
                continue
            
            # Get data
            try:
                data[node, ii, :] = [float(val) for val in line.split(delimiter)]
            # Must be column headings
            except ValueError:
                labels[node, :] = [label_prefix + heading for heading in line.split(delimiter)]
            
            ii += 1
    
    
    
    #%% Plot
    for node_idx in range(node_start, node_end):
        for col_idx in range(col_start, col_end):
            if data.shape[1] == 1:
                plt.plot(data[node_idx, :, 0], label=labels[node_idx][col_idx])
            else:
                plt.plot(data[node_idx, :, 0], data[node_idx, :, col_idx], label=labels[node_idx][col_idx])
    
    plt.xlabel(labels[0,0])
    plt.legend()
    plt.show()



#%% Run main()
if __name__ == "__main__":
    main(sys.argv[1], **dict(arg.split('=') for arg in sys.argv[2:]))