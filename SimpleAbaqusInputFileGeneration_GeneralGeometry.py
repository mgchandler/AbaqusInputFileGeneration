# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 09:34:36 2016

@author: ab9621
"""

import numpy as np
import os
import subprocess

def node_locs(nodes):
    #Converts the locations of nodes passed as strings into floats
    return (np.array([a.split() for a in nodes[1:-2]])).astype(float)
    
def find_nodes_in_x_width(nodes, target, width):
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
    # print(count)
    return locs[:count]   

def find_nearest_nodes(nodes, targets):
    #Find and return the nearest node to the target. Also
    #currently assumes the transducer is on the top surface.
    nodes = node_locs(nodes)
    nearest_nodes = np.ones((targets.shape[1]+1, len(nodes)))
    count = 0
    for a in nodes:
        if a[2] == 0.:
            nearest_nodes[0, count] = int(a[0])
            for el in range(targets.shape[1]):
                nearest_nodes[el+1, count] = np.abs(targets[0, el] - a[1])
            count += 1
    nearest_logical = np.argmin(nearest_nodes, 1)
    return nearest_nodes[0, nearest_logical[1:]].astype(int)
    
def write_input():
    path = 'C:\\Users\\mc16535\\OneDrive - University of Bristol\\Documents\\Postgrad\\Coding\\Abaqus\\FMC Generation\\v2'
    os.chdir(path)
    
    #set some global bits and bobs
    job_name_template = 'FMC_32els_rect_SDH'
    
    #Material Properties
    material_name = 'Aluminium'
    density = 2700.
    modulus = 70e9
    poisson = 0.34
    v = np.sqrt(modulus * (1 - poisson) / (density * (1 + poisson) * (1 - 2*poisson)))
    
    #Probe properties
    f = 5E6
    wavelength = v/f
    probe_width = 0.95e-3
    probe_pitch = probe_width + 0.05e-3
    num_els = 32
    
    amplitude = np.loadtxt('5MHz pulse.csv', delimiter=',')
    
    #Mesh properties
    n_per_wavelength = 20.
    dx = wavelength/n_per_wavelength
    print('dx = {}'.format(dx))
    element_type = 'CPS3'
    
    #Get probe centres
    probe_centres = np.zeros((2, num_els))
    probe_centres[0, :] = np.linspace(0, probe_pitch*(num_els-1), num_els)
    probe_centres[0, :] -= np.mean(probe_centres[0, :])
    
    #Generate a set of points to discretise transducer surfaces.
    num_nodes_on_tx = np.round(probe_width/dx).astype(int)
    probe_coords = np.zeros((2, num_els, num_nodes_on_tx))
    for el in range(num_els):
        probe_coords[0, el, :] = np.linspace(
            probe_centres[0, el]-probe_width/2,
            probe_centres[0, el]+probe_width/2,
            num_nodes_on_tx
        )
    probe_coords2 = np.reshape(probe_coords, (2, num_nodes_on_tx*num_els))
    
    #Work out the maximum allowed area for each triangle
    max_area = (dx*dx)/2.
    #max_area = 0.01
    print('max area = {}\n'.format(max_area))
    
    #Create the shape geometry. List coordinates adjacently starting with the
    #one on the +ve x-axis.
    ext_corners = np.array([
        [0.1,  0.1,  -0.1, -0.1],
        [0.,-20e-3,-20e-3,  0.]
    ])
    
    #Step properties
    time_step = 2 * dx/v
    step_time_length = 1.75 * np.max(np.abs(ext_corners[0, :]))/v
    
    ext_corners = np.append(probe_centres, ext_corners, 1)
    outer_vertices = ext_corners.shape[1]
    
    #Create scattering object
    xS = [ 27.5e-3]
    yS = [-15.0e-3]
    rS = [  1.5e-3]
    
    for hole in range(len(xS)):
        N = int(np.round(2*np.pi*rS[hole] / dx))
        jj = np.linspace(1, N, N)
        SDH = np.exp(2*np.pi*1j * jj/N) * rS[hole]
        x_full = np.append(ext_corners[0, :], np.real(SDH) + xS[hole])
        y_full = np.append(ext_corners[1, :], np.imag(SDH) + yS[hole])
    ext_corners = np.array([x_full, y_full])
        
    hole_vertices = SDH.shape[0]
    
    count = 1
    
    #create the poly file for input into Triangle
    filename = '{}.poly'.format(job_name_template)
    
    #Write out the node numbers and locations
    with open(filename, 'w') as f:
        f.write('{} 2 0 0\n'.format(ext_corners.shape[1])) #num points, 2 dimensions, no attributes, no boundary markers
        for ii in range(ext_corners.shape[1]):
            f.write('{}   {} {}\n'.format(ii+1, ext_corners[0, ii], ext_corners[1, ii]))# str(count)+ '   ' + str(a) + ' ' + str(b) +'\n')
        
        #Write out the number of segments and number of boundary markers
        f.write('{} 1\n'.format(ext_corners.shape[1])) #number of segments, number of boundary markers
        count = 1
        for ii in range(ext_corners.shape[1] - N):
            f.write('{}    {} {} {}\n'.format(count, (ii+1), (ii+1)%outer_vertices+1, 1))
            count += 1
        for ii in range(N):
            f.write('{}    {} {} {}\n'.format(count, ii+1 + ext_corners.shape[1] - N, (ii+1)%hole_vertices+1 + ext_corners.shape[1] - N, 2))
            count += 1
        f.write('{}\n'.format(len(xS)))
        if len(xS) >= 1:
            for hole in range(len(xS)):
                f.write('{}   {} {}\n'.format(hole+1, xS[hole], yS[hole]))
                
                
    
    #Call Triangle to generate the tri mesh
    subprocess.call("triangle -p -a{:.15f} -j -q30 -C {}".format(max_area, filename))
    
    f = open('{}.1.node'.format(job_name_template), 'r')
    nodes = f.read().split('\n')
    node_lines = len(nodes)
    n_nodes = node_lines-3
    f.close()
    f = open('{}.1.ele'.format(job_name_template), 'r')
    elements = f.read().split('\n')
    ele_lines = len(elements)
    n_elements = ele_lines-3
    f.close()
    
    print('{}'.format(job_name_template))
    
    for el in range(num_els):
        
        input_name = '{}_{}.inp'.format(job_name_template, el+1)
        
        #Open the input file and write the preamble
        with open(input_name, 'w') as i:
            i.write('*Heading\n')
            i.write('** Job name: {} Model name: Model-1\n'.format(job_name_template))
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
            
            #Create the all nodes and all elements sets
            i.write('*Nset, nset=All_nodes, generate\n')
            i.write('1, {}, 1\n'.format(n_nodes))
            
            i.write('*Elset, elset=All_elements, generate\n')
            i.write('1, {}, 1\n'.format(n_elements))
            
            #Create the node set which corresponds to the transmitting transducer
            i.write('*Nset, nset=FullProbe\n')
            for a in range(num_els):
                i.write('{},\n '.format(a+1))
            
            #Create the node set which corresponds to the transmitting transducer
            i.write('*Nset, nset=Transducer\n')
            transducer_nodes = find_nodes_in_x_width(nodes, probe_centres[:, el], probe_width)
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
            # i.write('*Output, field, time interval={:2.1g}\n'.format(step_time_length/100))
            # i.write('*Node Output\n')
            # i.write('POR,\n')
            # i.write('*Element Output\n')
            # i.write('S,\n')
            i.write('*Output, history, frequency=2\n')
            i.write('*Node Output, nset=FullProbe\n')
            i.write('U2,\n')
            
            #End of step
            i.write('*End Step\n')
            
            print('Tx = {} written'.format(el+1))
    
    print('Input file written')
    
if __name__ == '__main__':
    write_input()