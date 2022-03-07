# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 09:34:36 2016

@author: ab9621
"""

import numpy as np
import os
import subprocess
import time

def node_locs(nodes):
    #Converts the locations of nodes passed as strings into floats
    return (np.array([a.split() for a in nodes[1:-2]])).astype(float)
    
def find_nodes_at_loc(nodes, target, width):
    #Find and return the nodes at the target in a given width
    #Its to find the nodes in a transducer face for 2D
    #Currently assuming the transducer is on the top surface
    nodes = node_locs(nodes)
    locs = [0]*1000 #arbitrary long list
    count_1 = 0
    count_2 = 0
    if len(target.shape) == 1:
        x_min = target[0] - width/2.
        x_max = target[0] + width/2.
    else:
        target_xs = target[0, :]
        target = target[:, 0]
        nearest = np.ones((len(target_xs)+1, len(nodes)))
    for a in nodes:
        if a[2] == target[1]:
            if width != 0 and x_min <= a[1] <= x_max:
                locs[count_1] = int(a[0])
                count_1 += 1
            # In the event that there are no nodes within the width, find the
            # closest node instead.
            elif width == 0:
                nearest[0, count_2] = int(a[0])
                for el in range(len(target_xs)):
                    nearest[el+1, count_2] = np.abs(target_xs[el] - a[1])
                count_2 += 1
    if width != 0:
        return locs[:count_1] 
    else:
        nearest_logical = np.argmin(nearest, 1)
        return nearest[0, nearest_logical[1:]].astype(int)
    
def write_input():
    
    #set some global bits and bobs
#    path = 'O:/Documents/Work/Pogo FEA/pogo 1-1-164 learning/Writing input files using python'
    path = 'C:\\Users\\mc16535\\OneDrive - University of Bristol\\Documents\\Postgrad\\Coding\\Abaqus\\FMC Generation\\v1\\io files'
    os.chdir(path)
    
    #Material Properties
    material_name = 'Aluminium'
    v = 6400.
    density = 2700.
    modulus = 7e+10
    poisson = 0.34
    
    #Probe properties
    f = 5E6
    wavelength = v/f
    el_width = 1.95e-3
    el_sep = 0.05e-3
    el_pitch = el_width + el_sep
    probe_els = 16
    el_centres = np.zeros((2, probe_els))
    el_centres[0, :] = np.linspace(0, el_pitch*(probe_els-1), probe_els)
    el_centres[0, :] -= np.mean(el_centres[0, :])
    
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
    time_step = 1e-8
    step_time_length = 7.5e-5
    
    #Create the shape geometry
    x_outer = np.array([ -40e-3,  40e-3,  40e-3, -40e-3])
    y_outer = np.array([   0e-3,   0e-3, -40e-3, -40e-3])
    outer_vertices = len(x_outer)
    
    #Create scattering object
    xS = [ 27.5e-3]
    yS = [-30.0e-3]
    rS = [  1.5e-3]
    
    for hole in range(len(xS)):
        N = int(np.round(2*np.pi*rS[hole] / dx))
        jj = np.linspace(1, N, N)
        SDH = np.exp(2*np.pi*1j * jj/N) * rS[hole]
        x_full = np.append(x_outer, np.real(SDH) + xS[hole])
        y_full = np.append(y_outer, np.imag(SDH) + yS[hole])
        
    hole_vertices = len(SDH)
        
    
    
    #create the poly file for input into Triangle
    filename_template = '16El_FMC'
    filename = '{}.poly'.format(filename_template)
    
    #Write out the node numbers and locations
    f = open(filename, 'w')
    f.write('{} 2 0 0\n'.format(len(x_full))) #4 points, 2 dimensions, no attributes, no boundary markers
    for ii in range(len(x_full)):
        f.write('{}   {} {}\n'.format(ii+1, x_full[ii], y_full[ii]))# str(count)+ '   ' + str(a) + ' ' + str(b) +'\n')
    
    #Write out the number of segments and number of boundary markers
    f.write('{} 1\n'.format(len(x_full))) #number of segments, number of boundary markers
    count = 1
    for ii in range(len(x_full) - N):
        f.write('{}    {} {} {}\n'.format(count, (ii+1), (ii+1)%outer_vertices+1, 1))
        count += 1
    for ii in range(N):
        f.write('{}    {} {} {}\n'.format(count, (ii+1) + len(x_full) - N, (ii+1)%hole_vertices+1 + len(x_full) - N, 2))
        count += 1
    f.write('{}\n'.format(len(xS)))
    if len(xS) >= 1:
        for hole in range(len(xS)):
            f.write('{}   {} {}\n'.format(hole+1, xS[hole], yS[hole]))
    
    f.close()
    
    #Call Triangle to generate the tri mesh
    subprocess.call("triangle -p -a{:.15f} -j -q30 -C {}".format(max_area, filename))
    
    
    
    f = open('{}.1.node'.format(filename_template), 'r')
    nodes = f.read().split('\n')
    node_lines = len(nodes)
    n_nodes = node_lines-3
    f.close()
    f = open('{}.1.ele'.format(filename_template), 'r')
    elements = f.read().split('\n')
    ele_lines = len(elements)
    n_elements = ele_lines-3
    f.close()
    
    full_probe_nodes = find_nodes_at_loc(nodes, el_centres, 0)
        
    
    
    for element in range(probe_els):
        
        probe_centre = el_centres[:, element] #centre in global coordinates
        
        filename_el_template = '{}_el_{}'.format(filename_template, element)
    
        count = 1
        input_name = '{}.inp'.format(filename_el_template)
        job_name = '{}'.format(filename_el_template)
        
        #Open the input file and write the preamble
        i = open(input_name, 'w')
        i.write('*Heading\n')
        i.write('** Job name: {} Model name: Model-1\n'.format(job_name))
        i.write('*Preprint, echo=NO, model=NO, history=NO, contact=NO\n')
        i.write('**\n')
        i.write('** PART INSTANCE: Part-1-1\n')
        
        #Write the node locations
        i.write('*Node\n')
        
        for c1 in range(1,node_lines-2): #1 and -2 get rid of the extra bits
            line = list(filter(None, (nodes[c1]).split(' ')))
            i.write(', '.join(line[:3]) + '\n')
        
        #Write the element types
        i.write('*Element, type={}\n'.format(element_type))
        
        for c1 in range(1,ele_lines-2): #1 and -2 get rid of the extra bits
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
        transducer_nodes = find_nodes_at_loc(nodes, probe_centre, el_width)
        for a in transducer_nodes:
            i.write('{},\n '.format(a))
        
        #Create the node set which corresponds to the other array elements
        i.write('*Nset, nset=Full_Probe\n')
        for el in range(probe_els):
            i.write('{},\n '.format(full_probe_nodes[el]))
        
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
        i.write('*Step, name=Step-1, nlgeom=NO\n')
        i.write('*Dynamic, Explicit, direct user control\n')
        i.write('{}, {}\n'.format(time_step, step_time_length))
        
        #Write bulk viscosity, don't know if this is needed
        i.write('*Bulk Viscosity\n')
        i.write('0.0, 0.0\n')
        
        #Loads
        i.write('*Cload, amplitude=Amp-1\n')
        i.write('Transducer, 2, -10\n')
        
        #Outputs
        # i.write('*Output, field, time interval={:2.1g}\n'.format(step_time_length/5))
        # i.write('*Node Output\n')
        # i.write('POR,\n')
        # i.write('*Element Output\n')
        # i.write('S,\n')
        i.write('*Output, history, frequency=2\n')
        i.write('*Node Output, nset=Full_Probe\n')
        i.write('U2,\n')
        
        #End of step
        i.write('*End Step\n')
        
        i.close()
        
        print('El {} written'.format(element))
    
    print('Input files written')
    
if __name__ == '__main__':
    write_input()