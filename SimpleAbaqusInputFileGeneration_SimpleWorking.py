# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 09:34:36 2016

@author: ab9621
"""

import numpy as np
import subprocess

def node_locs(nodes):
    #Converts the locations of nodes passed as strings into floats
    return (np.array([a.split() for a in nodes[1:-2]])).astype(float)
    
def find_nodes_at_loc(nodes, target, width):
    #Find and return the nodes at the target in a given width
    #Its to find the nodes in a transducer face for 2D
    #Currently assuming the transducer is on the top surface
    nodes = node_locs(nodes)
    locs = [0]*1000 #arbitrary long list
    count = 0
    x_min = target[0] - width/2.
    x_max = target[0] + width/2.
    for a in nodes:
        if a[2] == target[1]:
            if x_min <= a[1] <= x_max:
                locs[count] = int(a[0])
                count += 1
    print(count)
    return locs[:count]                
    
def write_input():
    
    #set some global bits and bobs
#    path = 'O:/Documents/Work/Pogo FEA/pogo 1-1-164 learning/Writing input files using python'
    path = 'C:/Google Drive/pogo/input_files_python'
    input_name = 'Attempt1.inp'
    job_name = 'Attempt 1'
    
    #Material Properties
    material_name = 'Aluminium'
    v = 6400.
    density = 2700.
    modulus = 7e+10
    poisson = 0.34
    
    #Probe properties
    f = 5E6
    wavelength = v/f
    probe_centre = [0.025, 0.02] #centre in global coordinates
    probe_width = 5e-3
    
    amplitude = np.loadtxt('5MHz pulse.csv', delimiter=',')
    
    #Mesh properties
    n_per_wavelength = 20.
    dx = wavelength/n_per_wavelength
    print('dx = {}'.format(dx))
    element_type = 'CPS3'
    
    #Work out the maximum allowed area for each triangle
    max_area = (dx*dx)/2.
    #max_area = 0.01
    print('max area = {}\n'.format(max_area))
    
    #Step properties
    time_step = 3e-9
    step_time_length = 5e-5
    
    #Create the shape geometry
    x0 = 0.0
    x1 = 0.5
    
    y0 = 0.0
    y1 = 0.02
    
    x = [x0, x1]
    y = [y0, y1]
    count = 1
    
    #create the poly file for input into Triangle
    filename = 'Attempt1.poly'
    
    #Write out the node numbers and locations
    f = open(filename, 'w')
    f.write('4 2 0 0\n') #4 points, 2 dimensions, no attributes, no boundary markers
    for a in x:
        for b in y:
            f.write('{}   {} {}\n'.format(count, a, b))# str(count)+ '   ' + str(a) + ' ' + str(b) +'\n')
    
    #Write out the number of segments and number of boundary markers
    f.write('4 0\n') #number of segments, number of boundary markers
    count = 1
    f.write('1    1 2\n')
    f.write('2    2 4\n')
    f.write('3    4 3\n')
    f.write('4    3 1\n')
    f.write('0')
    
    f.close()
    
    #Open the input file and write the preamble
    i = open(input_name, 'w')
    i.write('*Heading\n')
    i.write('** Job name: {} Model name: Model-1\n'.format(job_name))
    i.write('*Preprint, echo=NO, model=NO, history=NO, contact=NO\n')
    i.write('**\n')
    i.write('** PART INSTANCE: Part-1-1\n')
    
    #Call Triangle to generate the tri mesh
    #subprocess.call("triangle -p -a{:.15f} -j -q30 -C {}".format(max_area, filename), cwd=path)
    subprocess.call("triangle -p -a{:.15f} -j -q30 -C {}".format(max_area, filename))
    #Write the node locations
    i.write('*Node\n')
    f = open('Attempt1.1.node', 'r')
    nodes = f.read().split('\n')
    n_lines = len(nodes)
    n_nodes = n_lines-3
    
    for c1 in range(1,n_lines-2): #1 and -2 get rid of the extra bits
        line = list(filter(None, (nodes[c1]).split(' ')))
        i.write(', '.join(line[:3]) + '\n')
    f.close()
    
    #Write the element types
    i.write('*Element, type={}\n'.format(element_type))
    f = open('Attempt1.1.ele', 'r')
    elements = f.read().split('\n')
    n_lines = len(elements)
    n_elements = n_lines-3
    
    for c1 in range(1,n_lines-2): #1 and -2 get rid of the extra bits
        line = filter(None, (elements[c1]).split(' '))
        i.write(', '.join(line) + '\n')    
    f.close()
    
    #Create the all nodes and all elements sets
    i.write('*Nset, nset=All_nodes, generate\n')
    i.write('1, {}, 1\n'.format(n_nodes))
    
    i.write('*Elset, elset=All_elements, generate\n')
    i.write('1, {}, 1\n'.format(n_elements))
    
    #Create the node set which corresponds to the transducer
    i.write('*Nset, nset=Transducer\n')
    transducer_nodes = find_nodes_at_loc(nodes, probe_centre, probe_width)
    for a in transducer_nodes:
        i.write('{},\n '.format(a))
    
    #Write the section definition
    i.write('*Solid Section, elset=All_elements, material={}\n'.format(material_name))
    i.write(',\n')
    
    #Write the amplitude definition
    i.write('*System\n')
    i.write('*Amplitude, name=Amp-1\n')
    #i.write('{}\n'.format(', '.join(map(str, amplitude))))
    for a in amplitude:
        i.write('{}\n'.format(', '.join(map(str, a))))
    
    #Write the material definition
    i.write('*Material, name={}\n'.format(material_name))
    i.write('*Density\n')
    i.write('{},\n'.format(density))
    i.write('*Elastic\n')
    i.write('{}, {}\n'.format(modulus, poisson))
    
    #Write the step
    i.write('*Step, name=Step-1, nlgeom=YES\n')
    i.write('*Dynamic, Explicit, direct user control\n')
    i.write('{}, {}\n'.format(time_step, step_time_length))
    
    #Write bulk viscosity, don't know if this is needed
    i.write('*Bulk Viscosity\n')
    i.write('0.0, 0.0\n')
    
    #Loads
    i.write('*Cload, amplitude=Amp-1\n')
    i.write('Transducer, 2, -10\n')
    
    #Outputs
    i.write('*Output, field, time interval={:2.1g}\n'.format(step_time_length/100))
    i.write('*Node Output\n')
    i.write('POR,\n')
    i.write('*Element Output\n')
    i.write('S,\n')
    i.write('*Output, history, frequency=2\n')
    i.write('*Node Output, nset=Transducer\n')
    i.write('U2,\n')
    
    #End of step
    i.write('*End Step\n')
    
    i.close()
    
    print('Input file written')
    
if __name__ == '__main__':
    write_input()