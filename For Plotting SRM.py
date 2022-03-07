# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 09:34:36 2016

@author: ab9621
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess
import sys
import yaml



def lame_const(E, nu):
    L = nu * E / ((1 + nu) * (1 - 2*nu))
    M = E / (2 * (1 + nu))
    return [L, M]

def velocity(E, nu, rho):
    v_L = np.sqrt(E * (1 - nu) / (rho * (1 + nu) * (1 - 2*nu)))
    v_S = np.sqrt(E / (2 * rho * (1 + nu)))
    return [v_L, v_S]

def wavelength(v, f):
    return v/f

def bigX(SRM_boundary, SRM_thickness, x, p):
    return ((SRM_boundary - x) / (SRM_thickness))**p

def attn(alpha_max, X):
    return alpha_max * X

def bigE(threshold, SRM_thickness, lam, SRM_boundary, E_0, x):
    k_i = 2*np.pi/lam
    alpha_max = - np.log(threshold) / (SRM_thickness * k_i)
    X = bigX(SRM_boundary, SRM_thickness, x, 3)
    return E_0 * np.exp(-alpha_max * X * k_i * SRM_thickness*bigX(SRM_boundary, SRM_thickness, x, 1))
    


def read_settings(filename):
    with open(filename, 'r') as f:
        settings = yaml.safe_load(f)
    return settings

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
    
def write_input(settings):
    
    #set some global bits and bobs
    job_name_template = settings['job']['name']
    
    #Material Properties
    material_name = settings['material']['name']
    density = settings['material']['density']
    modulus = settings['material']['modulus']
    poisson = settings['material']['poisson']
    v_L = velocity(modulus, poisson, density)[0]
    
    #Probe properties
    freq = settings['probe']['freq']
    lam_L = wavelength(v_L, freq)
    probe_width = settings['probe']['width']
    probe_pitch = probe_width + settings['probe']['separation']
    num_els = settings['probe']['num_els']
    
    amplitude = np.loadtxt(settings['probe']['amplitude'], delimiter=',')
    
    #Mesh properties
    n_per_wavelength = settings['mesh']['n_per_wl']
    dx = lam_L/n_per_wavelength
    print('dx = {}'.format(dx))
    element_type = settings['mesh']['el_type']
    
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
    probe_coords = np.reshape(probe_coords, (2, num_nodes_on_tx*num_els))
    
    #Work out the maximum allowed area for each triangle
    max_area = (dx*dx)/2.
    #max_area = 0.01
    print('max area = {}\n'.format(max_area))
    
    #Create the shape geometry. List coordinates adjacently starting with the
    #one on the +ve x-axis.
    ext_corners = np.array([
        settings['mesh']['geom']['x'],
        settings['mesh']['geom']['y']
    ])
    
    #Step properties
    time_step = 0.1 * dx/v_L
    step_time_length = settings['probe']['time_len']#1.75 * np.max(np.abs(ext_corners[0, :]))/v
    
    ext_corners = np.append(probe_coords, ext_corners, 1)
    
    
    
    #Create nodes for SRM layers. Segments will be drawn to separate layers, and
    # triangle will be allowed to discretise up each layer into elements.
    SRM_thickness = 1.5 * lam_L
    SRM_n_layers = int(np.round(SRM_thickness / (dx)) + 1)
    SRM_node_xs = np.linspace(ext_corners[0][-1], ext_corners[0][-1] - SRM_thickness, SRM_n_layers)
    SRM_nodes1 = np.array(
        [SRM_node_xs,
        np.full(SRM_node_xs.shape, ext_corners[1, -2])]
    )
    SRM_nodes2 = np.flip(np.array(
        [SRM_node_xs,
        np.full(SRM_node_xs.shape, ext_corners[1, -1])]
    ), axis=1)
    ext_corners = np.insert(ext_corners, [-1], SRM_nodes1[:, 1:], axis=1)
    ext_corners = np.insert(ext_corners, [-1], SRM_nodes2[:, :-1], axis=1)
    outer_vertices = ext_corners.shape[1]
    
    
    
    #Create scattering object
    xS = settings['mesh']['sdh']['x']
    yS = settings['mesh']['sdh']['y']
    rS = settings['mesh']['sdh']['r']
    
    for hole in range(len(xS)):
        N_hole = int(np.round(2*np.pi*rS[hole] / dx))
        jj = np.linspace(1, N_hole, N_hole)
        SDH = np.exp(2*np.pi*1j * jj/N_hole) * rS[hole]
        x_full = np.append(ext_corners[0, :], np.real(SDH) + xS[hole])
        y_full = np.append(ext_corners[1, :], np.imag(SDH) + yS[hole])
    ext_corners = np.array([x_full, y_full])
        
    hole_vertices = SDH.shape[0]
    
    #create the poly file for input into Triangle
    filename = '{}.poly'.format(job_name_template)
    node_segments = np.zeros((2, 10000)) #Arbitrary large number for plotting
    
    #Write out the node numbers and locations
    with open(filename, 'w') as f:
        f.write('{} 2 0 0\n'.format(ext_corners.shape[1])) #num points, 2 dimensions, no attributes, no boundary markers
        for ii in range(ext_corners.shape[1]):
            f.write('{}   {} {}\n'.format(ii+1, ext_corners[0, ii], ext_corners[1, ii]))# str(count)+ '   ' + str(a) + ' ' + str(b) +'\n')
        
        #Write out the number of segments and number of boundary markers
        f.write('{} 1\n'.format(ext_corners.shape[1]+SRM_n_layers-1)) #number of segments, number of boundary markers
        count = 1
        #Outer vertices
        for ii in range(ext_corners.shape[1] - N_hole):
            f.write('{}    {} {} {}\n'.format(count, (ii+1), (ii+1)%outer_vertices+1, 1))
            node_segments[:, count-1] = [(ii+1), (ii+1)%outer_vertices+1]
            count += 1
        #Hole vertices
        for ii in range(N_hole):
            f.write('{}    {} {} {}\n'.format(count, ii+1 + ext_corners.shape[1] - N_hole, (ii+1)%hole_vertices+1 + ext_corners.shape[1] - N_hole, 2))
            
            node_segments[:, count-1] = [ii+1 + ext_corners.shape[1] - N_hole, (ii+1)%hole_vertices+1 + ext_corners.shape[1] - N_hole]
            count += 1
        #SRM vertices
        SRM_start_idx = ext_corners.shape[1] - N_hole - 2*(SRM_n_layers - 1) + 1
        for ii in range(SRM_n_layers-1):
            f.write('{}    {} {} {}\n'.format(count, SRM_start_idx+ii-2, SRM_start_idx+2*SRM_n_layers-3-ii, 3)) #Maybe change these back to 3, 4 if meshing doesn't work in ABQ
            node_segments[:, count-1] = [SRM_start_idx+ii-2, SRM_start_idx+2*SRM_n_layers-3-ii]
            count += 1
        f.write('{}\n'.format(len(xS)))
        if len(xS) >= 1:
            for hole in range(len(xS)):
                f.write('{}   {} {}\n'.format(hole+1, xS[hole], yS[hole]))
                
                
                
    node_segments = node_segments[:, :389]
    
    min_x = np.min(ext_corners[0, :])
    max_x = np.max(ext_corners[0, :])
    min_y = np.min(ext_corners[1, :])
    max_y = np.max(ext_corners[1, :])
    Lx = max_x - min_x
    Ly = max_y - min_y
    L = np.max([Lx, Ly])
    num_E_pts = int(1.1 * L / (2*dx))
    
    E_xs = np.linspace(min_x-Lx*0.1, max_x+Lx*0.1, num_E_pts)
    E_ys = np.linspace(min_y-Ly*0.1, max_y+Ly*0.1, num_E_pts)
    E_vals = modulus*1e-9 * np.zeros((num_E_pts, num_E_pts))
    
    for jj in range(num_E_pts):
        yy = E_ys[jj]
        
        # If we are above or below the domain
        if (yy > 3e-16) or (yy + 45e-3 < 3e-16):
            E_vals[:, jj] = np.nan
        # Then we are inside the domain
        else:
            for ii in range(num_E_pts):
                xx = E_xs[ii]
                # If we are below the primary frontwall
                # (i.e. we are in the well on the RHS)
                if yy + 25e-3 < 3e-16:
                    if xx - 25e-3 < 3e-16:
                        E_vals[ii, jj] = np.nan
                    elif xx - 40e-3 > 3e-16:
                        E_vals[ii, jj] = np.nan
                # If we are above the primary frontwall
                else:
                    if xx - min_x < 3e-16:
                        E_vals[ii, jj] = np.nan
                    elif xx - 40e-3 > 3e-16:
                        E_vals[ii, jj] = np.nan
                    # Are we inside the SDH?
                    elif (xx - xS[0])**2 + (yy - yS[0])**2 < rS[0]**2:
                        E_vals[ii, jj] = np.nan
                    # Are we in the SRM?
                    elif xx + 30e-3 < 3e-16:
                        # E_vals[ii, jj] = bigE(0.01, abs(min_x + 30e-3), lam_L, -30e-3, modulus, xx) * 1e-9
                        E_vals[ii, jj] = 2*np.pi*freq * bigX(-30e-3, abs(min_x + 30e-3), xx, 3)
                    
                
    
    fs = 12
    plt.figure(figsize=(8, 5), dpi=750)
    plt.imshow(np.transpose(E_vals), extent=[E_xs[0]*1e3, E_xs[-1]*1e3, E_ys[0]*1e3, E_ys[-1]*1e3], origin='lower', cmap='viridis', aspect=.1)
    plt.colorbar(label=r"Rayleigh Mass-Damping Coefficient $C_M$")
    # plt.scatter(ext_corners[0, :], ext_corners[1, :], s=0.2)
    for ii in range(node_segments.shape[1]):
        if node_segments[0][ii] != 0 and node_segments[0][ii] != 0:
            plt.plot([ext_corners[0][int(node_segments[0][ii])-1]*1e3, ext_corners[0][int(node_segments[1][ii])-1]*1e3], 
                      [ext_corners[1][int(node_segments[0][ii])-1]*1e3, ext_corners[1][int(node_segments[1][ii])-1]*1e3],
                      'r', linewidth=0.75)
    plt.xlim(-33, -29)
    plt.ylim(-27, 2)
    plt.xlabel('x (mm)')
    plt.ylabel('z (mm)')
    plt.show()
        
                
    
    a = 1
    
    # #Call Triangle to generate the tri mesh.
    # # If triangle keeps bringing up errors, check the number of nodes and the
    # # number of segments being reported in filename.poly
    # if sys.platform == 'win32':
    #     print("triangle -p -a{:.15f} -j -q30 -C {}".format(max_area, filename))
    #     subprocess.run("triangle -p -a{:.15f} -j -q30 -C {}".format(max_area, filename))
    # elif sys.platform == 'linux':
    #     os.system("./triangle -p -a{:.15f} -j -q30 -C {}".format(max_area, filename))
    
    # # Get all node locations
    # f = open('{}.1.node'.format(job_name_template), 'r')
    # nodes = f.read().split('\n')
    # node_lines = len(nodes)
    # n_nodes = node_lines-3
    # f.close()
    # # Get all element locations
    # f = open('{}.1.ele'.format(job_name_template), 'r')
    # elements = f.read().split('\n')
    # ele_lines = len(elements)
    # n_elements = ele_lines-3
    # f.close()
    
    # # Store nodes in an easy to read format for PML area calculation.
    # node_block = np.zeros((2, n_nodes))
    # for node in range(n_nodes):
    #     line = nodes[node+1].split('  ')
    #     line = [item for item in line if item != '']
    #     node_block[:, node] = [float(line[1]), float(line[2])]
        
    # # Store elements in an easy to read format
    # element_block = np.zeros((3, n_elements))
    # SRM_elements = np.zeros((SRM_n_layers-1, 10000), dtype=int) # Arbitrary large number to store elements.
    # All_but_SRM_els = np.zeros(n_elements, dtype=int)
    # SRM_els_ticker = np.zeros(SRM_n_layers-1, dtype=int)
    # Not_ticker = 0
    # count = 0
    # for element in range(n_elements):
    #     line = elements[element+1].split('  ')
    #     line = [item for item in line if item != '']
    #     element_block[:, element] = [int(line[1]), int(line[2]), int(line[3])]
        
    #     SRM_ticker = 0
    #     el_pos = 0
        
    #     #Check whether this element is in the SRM.
    #     for node in range(3):
    #         # Nodes are recorded by Triangle from index 1, python refers from index 0.
    #         # Contents of element_block are recorded in Triangle notation.
    #         node_pos = node_block[:, int(element_block[node, element])-1]
    #         # Want this to be <= 0. Use machine epsilon in case of equality.
    #         if node_pos[0] - SRM_node_xs[0] < 3e-16:
    #             SRM_ticker += 1
    #             el_pos += node_pos[0] / 3
    #         # If any nodes are not, then this element is not in the SRM. Skip the rest
    #         # of the nodes in this element.
    #         else:
    #             break
            
    #     # If all nodes are in the SRM, check which sublayer it is part of
    #     if SRM_ticker == 3:
    #         # The element must be at least in the 0th sublayer.
    #         which_sublayer = 0
    #         # Check if it is in the next one.
    #         while el_pos < SRM_node_xs[which_sublayer+1]:
    #             # If it is, update the current sublayer.
    #             which_sublayer += 1
            
    #         # We now know which sublayer this element is in. Add it to that
    #         # sublayer
    #         SRM_elements[which_sublayer, SRM_els_ticker[which_sublayer]] = int(line[0])
    #         SRM_els_ticker[which_sublayer] += 1
    #     else:
    #         All_but_SRM_els[Not_ticker] = int(line[0])
    #         Not_ticker += 1
            
    # SRM_elements = SRM_elements[:, :np.max(SRM_els_ticker)] #Note there are still some zeros here
    # All_but_SRM_els = All_but_SRM_els[All_but_SRM_els != 0]
    # # All_but_SRM_els = np.delete(All_but_SRM_els, (All_but_SRM_els==0))
    
            
            
            
    # # Get Young's modulus and Rayleigh Mass Damping values in SRM
    # SRM_boundary = settings['mesh']['geom']['x'][-1]
    # SRM_modulus = bigE(0.01, SRM_thickness, lam_L, SRM_boundary, modulus, SRM_node_xs)
    # SRM_CM = 2*np.pi*freq * bigX(SRM_boundary, SRM_thickness, SRM_node_xs, 3)
    
    
    
    # print('{}'.format(job_name_template))
    
    # Approx_time_per_el_per_step = 1.03e-7
    # Est_runtime = n_elements * step_time_length * Approx_time_per_el_per_step / time_step
    # print('Est runtime per transmitter: {:4.3g} s'.format(Est_runtime))
    
    # for el in range(num_els):
        
    #     input_name = '{}_{}.inp'.format(job_name_template, el+1)
        
    #     #Open the input file and write the preamble
    #     with open(input_name, 'w') as i:
    #         i.write('*Heading\n')
    #         i.write('** Job name: {} Model name: Model-1\n'.format(job_name_template))
    #         i.write('*Preprint, echo=NO, model=NO, history=NO, contact=NO\n')
    #         i.write('**\n')
    #         i.write('** PART INSTANCE: Part-1-1\n')
            
    #         #Write the node locations
    #         i.write('*Node\n')
    #         for c1 in range(1,node_lines-2): #1 and -2 get rid of the extra bits
    #             line = list(filter(None, (nodes[c1]).split(' ')))
    #             i.write(', '.join(line[:3]) + '\n')
            
    #         #Write the element types
    #         i.write('*Element, type={}\n'.format(element_type))
    #         for c1 in range(1,ele_lines-2): #1 and -2 get rid of the extra bits
    #             line = filter(None, (elements[c1]).split(' '))
    #             i.write(', '.join(line) + '\n')    
            
    #         #Create the all nodes and all elements sets
    #         i.write('*Nset, nset=All_nodes, generate\n')
    #         i.write('1, {}, 1\n'.format(n_nodes))
            
    #         i.write('*Elset, elset=All_elements, generate\n')
    #         i.write('1, {}, 1\n'.format(n_elements))
            
    #         #Create the node set which corresponds to the transmitting transducer
    #         i.write('*Nset, nset=FullProbe\n')
    #         for a in range(num_els*num_nodes_on_tx):
    #             i.write('{},\n '.format(a+1))
            
    #         #Create the node set which corresponds to the transmitting transducer
    #         i.write('*Nset, nset=Transducer\n')
    #         for a in range(num_nodes_on_tx):
    #             b = el * num_nodes_on_tx + a +1
    #             i.write('{},\n '.format(b))
                
    #         # Create element set for the SRM
    #         i.write('*Elset, elset=All_but_SRM\n')
    #         for a in All_but_SRM_els:
    #             i.write('{},\n'.format(a))
    #         for sublayer in range(1, SRM_n_layers):
    #             i.write('*Elset, elset=SRM{}\n'.format(sublayer))
    #             for a in SRM_elements[sublayer-1]:
    #                 if a != 0:
    #                     i.write('{}, \n'.format(a))
            
    #         #Write the section definition
    #         i.write('*Solid Section, elset=All_but_SRM, material={}\n'.format(material_name))
    #         i.write(',\n')
    #         for sublayer in range(1, SRM_n_layers):
    #             i.write('*Solid Section, elset=SRM{}, material={}{}\n'.format(sublayer, material_name, sublayer))
    #             i.write(',\n')
            
    #         #Write the amplitude definition
    #         i.write('*System\n')
    #         i.write('*Amplitude, name=Amp-1\n')
    #         #i.write('{}\n'.format(', '.join(map(str, amplitude))))
    #         for a in amplitude:
    #             i.write('{}\n'.format(', '.join(map(str, a))))
            
    #         #Write the material definition
    #         i.write('*Material, name={}\n'.format(material_name))
    #         i.write('*Density\n')
    #         i.write('{},\n'.format(density))
    #         i.write('*Elastic\n')
    #         i.write('{}, {}\n'.format(modulus, poisson))
    #         for sublayer in range(1, SRM_n_layers):
    #             i.write('*Material, name={}{}\n'.format(material_name, sublayer))
    #             i.write('*Density\n')
    #             i.write('{},\n'.format(density))
    #             i.write('*Elastic\n')
    #             i.write('{}, {}\n'.format(SRM_modulus[sublayer], poisson))
    #             i.write('*Damping, alpha={}\n'.format(SRM_CM[sublayer]))
                
            
    #         #Write the PML
    #         # if 'PML' in settings['mesh']:
    #         #     i.write('*PERFECTLY MATCHED LAYER, elset=PML, name=PMLayer, type=CARTESIAN\n')
    #         #     i.write('{}, {}, {}\n'.format(ext_corners[0, PML_idxs[0]-1], 0, ext_corners[1, PML_idxs[0]-1]))
    #         #     i.write('{}, {}, {}\n'.format(ext_corners[0, PML_idxs[1]-1], 0, ext_corners[1, PML_idxs[1]-1]))
    #             # i.write('*PML COEFFICIENT, variation=LINEAR\n')
            
    #         #Write the step
    #         i.write('*Step, name=Step-1, nlgeom=NO\n')
    #         i.write('*Dynamic, Explicit, direct user control\n')
    #         i.write('{}, {}\n'.format(time_step, step_time_length))
            
    #         #Write bulk viscosity, don't know if this is needed
    #         i.write('*Bulk Viscosity\n')
    #         i.write('0.0, 0.0\n')
            
    #         #Loads
    #         i.write('*Cload, amplitude=Amp-1\n')
    #         i.write('Transducer, 2, -10\n')
            
    #         #Outputs
    #         if 'field' in settings['output'].keys():
    #             i.write('*Output, field, time interval={}\n'.format(settings['output']['field']['t_int']))
    #             i.write('*{} Output, elset={}\n'.format(settings['output']['field']['output'], settings['output']['field']['elset']))
    #             i.write('{},\n'.format(settings['output']['field']['type']))
    #         if 'history' in settings['output'].keys():
    #             i.write('*Output, history, frequency={}\n'.format(settings['output']['history']['freq']))
    #             i.write('*{} Output, nset={}\n'.format(settings['output']['history']['output'], settings['output']['history']['nset']))
    #             i.write('{},\n'.format(settings['output']['history']['type']))
    #         # i.write('*Output, field, time interval={:2.1g}\n'.format(step_time_length/100))
    #         # i.write('*Node Output\n')
    #         # i.write('POR,\n')
    #         # i.write('*Element Output\n')
    #         # i.write('S,\n')
    #         # i.write('*Output, history, frequency=2\n')
    #         # i.write('*Node Output, nset=FullProbe\n')
    #         # i.write('U2,\n')
            
    #         #End of step
    #         i.write('*End Step\n')
            
    #         print('Tx = {} written'.format(el+1))
    
    # print('Input file written')
    
if __name__ == '__main__':
    if len(sys.argv) > 1:
        yaml_name = '{}.yml'.format(sys.argv[1])
    else:
        yaml_name = 'L_Vis.yml'
        os.chdir('C:\\Users\\mc16535\\OneDrive - University of Bristol\\Documents\\Postgrad\\Coding\\Abaqus\\FMC Generation\\v5')
    settings = read_settings(yaml_name)
    write_input(settings)