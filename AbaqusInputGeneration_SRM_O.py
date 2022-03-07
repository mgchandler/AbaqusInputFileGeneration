# -*- coding: utf-8 -*-
"""
Adapted on Tue Jun 18 17:28:19 2021

@author: mc16535

Original script created on Wed Jan 20 09:34:36 2016

#@author: ab9621
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess
import sys
import time
import yaml



#%% Support functions for wave and SRM parameters.

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

def bigX(SRM_n_layers, p):
    return np.linspace(0, 1, SRM_n_layers)**p

def attn(alpha_max, X):
    return alpha_max * X

def bigE(threshold, SRM_thickness, lam, SRM_n_layers, E_0):
    k_i = 2*np.pi/lam
    alpha_max = - np.log(threshold) / (SRM_thickness * k_i)
    X = bigX(SRM_n_layers, 3)
    return E_0 * np.exp(-alpha_max * X * k_i * SRM_thickness*bigX(SRM_n_layers, 1))
    


#%% Support functions for writing and node locations

def read_settings(filename):
    with open(filename, 'r') as f:
        settings = yaml.safe_load(f)
    return settings

def node_locs(nodes):
    # Converts the locations of nodes passed as strings into floats
    return (np.array([a.split() for a in nodes[1:-2]])).astype(float)
    
def find_nodes_in_x_width(nodes, target, width):
    # Find and return the node indices located within some horizontal width of
    # the target coordinates. This function is built on the assumption that the
    # transducer is oriented horizontally (parallel to x-axis), i.e. on a flat
    # top surface of a sample.
    
    nodes = node_locs(nodes)
    locs = [] # Initialise empty array: total # nodes in transducer unknown.
    x_min = target[0] - width/2.
    x_max = target[0] + width/2.
    for a in nodes:
        if a[2] == target[1]:
            if x_min <= a[1] <= x_max:
                locs.append(int(a[0]))
    return locs 

def checkIfLeft(linePoint1, linePoint2, checkPoint):
    return ((linePoint2[0] - linePoint1[0])*(checkPoint[1] - linePoint1[1]) - (linePoint2[1] - linePoint1[1])*(checkPoint[0] - linePoint1[0])) > -3e-16



#%% Generate nodes and write inputs
    
def write_input(settings):
    
    #%% Set global parameters
    
    job_name_template = settings['job']['name']
    
    # Radii and filenames
    meas_R1 = 10.0e-3
    meas_R2 = 20.0e-3
    meas_R3 = 30.0e-3
    rad_at_SRM = 32.5e-3
    topwallhalflen = 19e-3
    use_topwall = False
    
    # Material properties
    material_name = settings['material']['name']
    density = settings['material']['density']
    modulus = settings['material']['modulus']
    poisson = settings['material']['poisson']
    [v_L, v_S] = velocity(modulus, poisson, density)
    
    # Probe properties
    freq = settings['probe']['freq']
    lam_L = wavelength(v_L, freq)
    probe_width = settings['probe']['width']
    probe_pitch = probe_width + settings['probe']['separation']
    num_els = settings['probe']['num_els']
    
    amplitude = np.loadtxt(settings['probe']['amplitude'], delimiter=',')
    
    # Mesh properties
    n_per_wavelength = settings['mesh']['n_per_wl']
    dx = lam_L/n_per_wavelength # Approximate distance between nodes
    print('dx = {}'.format(dx))
    element_type = settings['mesh']['el_type']
    ext_corners = np.array([
        settings['mesh']['geom']['x'],
        settings['mesh']['geom']['y']
    ])
    
    geomR = rad_at_SRM + 1.5*lam_L
    N_wall = int(np.round(2*np.pi*geomR / dx))
    jj = np.linspace(1, int(N_wall), int(N_wall))
    
    outerDwall = np.exp(-2*np.pi*1j * jj/N_wall) * geomR
    # outerDwall[0] = geomR
    # outerDwall[-1] = -geomR
    ext_corners = np.array([
        np.real(outerDwall),
        np.imag(outerDwall)
    ])
    
    # Time step properties
    time_step = 0.1 * dx/v_L
    if 'time_len' in settings['probe'].keys():
        step_time_length = settings['probe']['time_len']
    else:        
        step_time_length = 1.1 * geomR / v_S
    
    # Get probe centre coordinates
    probe_centres = np.zeros((2, num_els))
    probe_centres[0, :] = np.linspace(0, probe_pitch*(num_els-1), num_els)
    probe_centres[0, :] -= np.mean(probe_centres[0, :])
    
    
    
    #%% Generate meshing algorithm inputs.
    # Node locations need to be determined.
    
    # Generate locations of nodes which will act as the transducer interfaces.
    # Make them regularly spaced so that the same loading amplitude can be
    # applied to them all.
    if probe_width == 0:
        num_nodes_on_tx = 1
    else:
        num_nodes_on_tx = np.round(probe_width/(.5*dx)).astype(int)
    probe_coords = np.zeros((2, num_els, num_nodes_on_tx))
    for el in range(num_els):
        probe_coords[0, el, :] = np.linspace(
            probe_centres[0, el]-probe_width/2,
            probe_centres[0, el]+probe_width/2,
            num_nodes_on_tx
        )
    N_probe_coords = num_nodes_on_tx*num_els
    probe_coords = np.reshape(probe_coords, (2, N_probe_coords))
    
    # Work out the maximum allowed area for each triangle
    max_area = (dx*dx)/2.
    print('max area = {}\n'.format(max_area))
    
    N_topwall = np.round(2*topwallhalflen/dx).astype(int)
    top_wall = np.zeros((2, N_topwall))
    top_wall[0, :] = np.linspace(-topwallhalflen, topwallhalflen, N_topwall)
    
    # Read the shape of the part to be inspected. List coordinates as they are
    # reached as the perimeter is traced out, starting in the +ve x-direction
    # from the probe.
    # if use_topwall:
    #     ext_corners = np.append(top_wall[:, 1:-1], ext_corners, 1)
    #     N_topwall -= 2
        
    #     N_surf = N_topwall
        
    #     assert(N_surf == amplitude.shape[1]-1)
    # else:
    #     ext_corners = np.append(probe_coords, ext_corners, 1)
        
    #     N_surf = N_probe_coords
    N_surf = 0
    
    
    
    # # Smooth out the corner from which an ultrasonic wave is diffracted. In the
    # # L-shaped geometry, this is usually the inner corner of the L, the 3rd index
    # # in the list of coordinates supplied.
    # if 'round' in settings['mesh'].keys():
    #     round_idx = settings['mesh']['round']['idx'] + probe_coords.shape[1]
    #     round_r = settings['mesh']['round']['r']
    #     N_round = int(np.round(2*np.pi*round_r / dx))
        
    #     if N_round > 0:
    #         # Get vectors of incoming wall and outgoing wall to fillet wall.
    #         vec1 = ext_corners[:, round_idx] - ext_corners[:, round_idx-1]
    #         vec1 = vec1 / np.linalg.norm(vec1)
    #         vec2 = ext_corners[:, round_idx+1] - ext_corners[:, round_idx]
    #         vec2 = vec2 / np.linalg.norm(vec2)
    #         # Angle incoming wall makes with the vertical.
    #         theta1 = np.arccos(np.dot([0, 1], vec1))
    #         # Angle the outgoing wall makes with the vertical.
    #         theta2 = np.arccos(np.dot([0, 1], vec2))
    #         # Half of the angle between incoming and outgoing wall.
    #         phi = np.arccos(np.dot(vec1, vec2))/2
    #         # Angle between the vertical and half way between incoming and outgoing
    #         # wall.
    #         alpha = theta1 + phi
            
    #         # Indices corresponding to the valid section of the circle to generate
    #         # the fillet.
    #         jj = np.linspace(1, N_round, N_round)
    #         jj = jj[int(np.round(theta1 * N_round/(2*np.pi))):int(np.round(theta2 * N_round/(2*np.pi)))]
    #         # Coordinates of the fillet, located in the right place by shifting
    #         # by the relative corner, and the radius of the curvature.
    #         round_ = np.exp(2*np.pi*1j * jj/N_round) * round_r + ext_corners[0, round_idx]+1j*ext_corners[1, round_idx] - round_r*np.sin(alpha)/np.sin(phi) - 1j*round_r*np.cos(alpha)/np.sin(phi) 
    #         round_ = np.array([np.real(round_), np.imag(round_)])
    #         N_round = round_.shape[1]-1
    #     else:
    #         round_ = np.reshape(ext_corners[:, round_idx], (2,1))
    #         N_round = round_.shape[1]-1
    #     ext_corners = np.concatenate((ext_corners[:, :round_idx], round_, ext_corners[:, round_idx+1:]), axis=1)
        
    
    
    # Stiffness Reduction Method layer (SRM layer; SRM; Pettit et al. 2014) is
    # used to absorb ultrasonic waves which are incident on walls that are
    # desired to be infinite in extent. The layer consists of several sub-layers
    # with decreasing Young's modulus and increasing damping as depth into the
    # SRM increases.
    # Pettit et al. 2014 DOI: 10.1016/j.ultras.2013.11.013
    # SRM was used as PML can only be used in Abaqus/Simple and not Abaqus/Explicit,
    # which has much longer time taken to compute per step. By using an SRM, we
    # increase setup time taken to run this .inp generation script (~5min -> 15
    # min), but massively reduce the time taken to run within Abaqus (~11hrs ->
    # 4hrs).
    # SRM in settings['mesh']['geom'] should list the indices of the corners which
    # denote SRM layers. Each element in the list n should be 0 <= n <= N-1 where
    # N is the total number of corners listed, and the boundary will be drawn along
    # the wall between corners n and n+1.
    N_SRM = 1
    SRM_thickness = 1.5 * lam_L
    SRM_n_layers = int(np.round(SRM_thickness / (dx)) + 1)
        
    SRM_nodes = np.array([np.linspace(rad_at_SRM, geomR, SRM_n_layers),
                          np.zeros((SRM_n_layers))])
    # ext_corners = np.insert(ext_corners, [N_surf], SRM_nodes[:, :-1], axis=1)
    # ext_corners = np.append(ext_corners, -np.flip(SRM_nodes[:, :-1], axis=1), axis=1)
    outer_vertices = ext_corners.shape[1]
    
    SRM_internal_nodes = []
    
    for layer in range(SRM_n_layers-1):
        SRM_R = SRM_nodes[0][layer]
        N_SRM = int(np.round(2*np.pi*SRM_R / dx))
        jj = np.linspace(1, int(N_SRM), int(N_SRM))
    
        innerDwalls = np.exp(-2*np.pi*1j * jj/N_SRM) * SRM_R
        innerDwalls = innerDwalls[1:-2]
        SRM_internal_nodes.append( np.array([
            np.real(innerDwalls),
            np.imag(innerDwalls)
        ]) )
    
    # if 'SRM' in settings['mesh']['geom']:
            
        
        
    #     boundaries = np.array(settings['mesh']['geom']['SRM'])
    #     b_counter = 0
    #     for boundary in boundaries:
    #         if 'round' not in settings['mesh'].keys():
    #             round_idx = 0
    #             N_round = 0
    #         assert boundary != round_idx - N_probe_coords
    #         if boundary < round_idx - N_probe_coords:
    #             boundaries[b_counter] = N_probe_coords + N_SRM + boundary
    #         elif boundary > round_idx - N_probe_coords:
    #             boundaries[b_counter] = N_probe_coords + N_round + 2*N_SRM + boundary
    #         boundary = boundaries[b_counter]
    #         b_counter += 1
            
    #         SRM_thickness = 1.5 * lam_L
    #         SRM_n_layers = int(np.round(SRM_thickness / (dx)) + 1)
    #         SRM_wall_dir = np.squeeze(ext_corners[:, boundary+1] - ext_corners[:, boundary])
    #         wall_theta = np.arccos(SRM_wall_dir[1] / np.sqrt(SRM_wall_dir[0]**2 + SRM_wall_dir[1]**2))
    #         SRM_nodes_out = np.array([np.linspace(ext_corners[0][boundary], ext_corners[0][boundary] - SRM_thickness*np.cos(wall_theta), SRM_n_layers),
    #                                   np.linspace(ext_corners[1][boundary], ext_corners[1][boundary] - SRM_thickness*np.sin(wall_theta), SRM_n_layers)])
    #         SRM_nodes_in =  np.array([np.linspace(ext_corners[0][boundary+1], ext_corners[0][boundary+1] - SRM_thickness*np.cos(wall_theta), SRM_n_layers),
    #                                   np.linspace(ext_corners[1][boundary+1], ext_corners[1][boundary+1] - SRM_thickness*np.sin(wall_theta), SRM_n_layers)])
    #         SRM_nodes = np.append(SRM_nodes_out, np.flip(SRM_nodes_in, axis=1), axis=1)
            
    #         N_SRM += int(SRM_nodes.shape[1]/2)-1
            
    #         # SRM layer is expected to be on the wall adjacent to the transducer surface
    #         # in the -x-direction, i.e. the penultimate wall in ext_corners.
    #         ext_corners = np.insert(ext_corners, [boundary+1], SRM_nodes[:, 1:-1], axis=1)
    #         outer_vertices = ext_corners.shape[1]
    
    # if 'sdh' in settings['mesh'].keys():
    #     # Create scattering object
    #     xS = settings['mesh']['sdh']['x']
    #     yS = settings['mesh']['sdh']['y']
    #     rS = settings['mesh']['sdh']['r']
        
    #     # Generate nodes for the boundary of the scatterer.
    #     for hole in range(len(xS)):
    #         N_hole = int(np.round(2*np.pi*rS[hole] / dx))
    #         jj = np.linspace(1, N_hole, N_hole)
    #         SDH = np.exp(2*np.pi*1j * jj/N_hole) * rS[hole]
    #         x_full = np.append(ext_corners[0, :], np.real(SDH) + xS[hole])
    #         y_full = np.append(ext_corners[1, :], np.imag(SDH) + yS[hole])
    #     ext_corners = np.array([x_full, y_full])
            
    #     hole_vertices = SDH.shape[0]
    # else:
    N_hole = 0
    
    N_meas1 = int(np.round(2*np.pi*meas_R1 / dx))
    jj = np.linspace(1, int(N_meas1), int(N_meas1))
    
    innerDmeas = np.exp(-2*np.pi*1j * jj/N_meas1) * meas_R1
    innerDmeas = innerDmeas[1:-2]
    meas_internal_nodes = np.array([
        np.real(innerDmeas),
        np.imag(innerDmeas)
    ])
    
    N_meas2 = int(np.round(2*np.pi*meas_R2 / dx))
    jj = np.linspace(1, int(N_meas2), int(N_meas2))
    
    innerDmeas = np.exp(-2*np.pi*1j * jj/N_meas2) * meas_R2
    innerDmeas = innerDmeas[1:-2]
    meas_internal_nodes = np.append(meas_internal_nodes, np.array([
        np.real(innerDmeas),
        np.imag(innerDmeas)
    ]), axis=1)
    
    N_meas3 = int(np.round(2*np.pi*meas_R3 / dx))
    jj = np.linspace(1, int(N_meas3), int(N_meas3))
    
    innerDmeas = np.exp(-2*np.pi*1j * jj/N_meas3) * meas_R3
    innerDmeas = innerDmeas[1:-2]
    meas_internal_nodes = np.append(meas_internal_nodes, np.array([
        np.real(innerDmeas),
        np.imag(innerDmeas)
    ]), axis=1)
    
    N_meas = N_meas1 + N_meas2 + N_meas3
    
    
    
    
    
    count = 1
    
    # .poly file used as input for Triangle (meshing algorithm). 
    # http://www.cs.cmu.edu/~quake/triangle.poly.html
    
    # Make two of them, one with the hole present and one without.
    job_name_template = '{}'.format(settings['job']['name'])
    filename = '{}.poly'.format(job_name_template)
    
    numNodes = ext_corners.shape[1]
    numNodes += probe_coords.shape[1]
    for layer in SRM_internal_nodes:
        numNodes += layer.shape[1]
    numNodes += meas_internal_nodes.shape[1]
    
    all_the_nodes = np.zeros((2, int(1.25*numNodes))) #Arbitrary large array; used for plotting only.
    node_segments = np.zeros((2, int(1.25*numNodes))) #Arbitrary large array; used for plotting only.
    
    # Write out the node numbers and locations
    with open(filename, 'w') as f:
        f.write('{} 2 0 0\n'.format(numNodes)) #num points, 2 dimensions, no attributes, no boundary markers
        # Write external corners
        count = 1
        for ii in range(ext_corners.shape[1]):
            f.write('{}   {} {}\n'.format(count, ext_corners[0, ii], ext_corners[1, ii]))
            all_the_nodes[:, count-1] = [ext_corners[0, ii], ext_corners[1, ii]]
            count += 1
        # Write probe coords
        for pr in range(probe_coords.shape[1]):
            f.write('{}   {} {}\n'.format(count, probe_coords[0, pr], probe_coords[1, pr]))
            all_the_nodes[:, count-1] = [probe_coords[0, pr], probe_coords[1, pr]]
            count += 1
        # Write SRM layers
        for layer in SRM_internal_nodes:
            for jj in range(layer.shape[1]):
                f.write('{}   {} {}\n'.format(count, layer[0, jj], layer[1, jj]))
                all_the_nodes[:, count-1] = [layer[0, jj], layer[1, jj]]
                count += 1
        # Write measuring layer
        for kk in range(meas_internal_nodes.shape[1]):
            f.write('{}   {} {}\n'.format(count, meas_internal_nodes[0, kk], meas_internal_nodes[1, kk]))
            all_the_nodes[:, count-1] = [meas_internal_nodes[0, kk], meas_internal_nodes[1, kk]]
            count += 1
        
        
        
        #Write out the number of segments and number of boundary markers
        numSegments = ext_corners.shape[1]
        for layer in SRM_internal_nodes:
            numSegments += layer.shape[1]
        
        f.write('{} 1\n'.format(numSegments)) #number of segments, number of boundary markers
        count = 1
        
        #Outer vertices (i.e. walls of geometry)
        for ii in range(ext_corners.shape[1] - N_hole):
            f.write('{}    {} {} {}\n'.format(count, (ii+1), (ii+1)%outer_vertices+1, 1))
            node_segments[:, count-1] = [(ii+1), (ii+1)%outer_vertices+1]
            count += 1
            
        #SRM vertices (pre-defined element edges. Not walls)
        N_uptolayer = ext_corners.shape[1] + probe_coords.shape[1]
        for ii in range(SRM_n_layers-1):
            layer = SRM_internal_nodes[ii]
            #Outer vertices (i.e. walls of geometry)
            for ii in range(layer.shape[1]):
                f.write('{}    {} {} {}\n'.format(count, N_uptolayer+(ii+1), N_uptolayer+(ii+1)%layer.shape[1]+1, 1))
                node_segments[:, count-1] = [N_uptolayer+(ii+1), N_uptolayer+(ii+1)%layer.shape[1]+1]
                count += 1
            N_uptolayer += layer.shape[1]
            # f.write('{}    {} {} {}\n'.format(count, N_surf+ii+1, count-ii, 3))
            # node_segments[:, count-1] = [N_surf+ii+1, count-ii]
            # count += 1
            # for jj in range(layer.shape[1]-1):
            #     f.write('{}    {} {} {}\n'.format(count, count-1-ii, count-ii, 3))
            #     node_segments[:, count-1] = [count-1-ii, count-ii]
            #     count += 1
            # f.write('{}    {} {} {}\n'.format(count, count-1-ii, ext_corners.shape[1]-ii, 3))
            # node_segments[:, count-1] = [count-1-ii, ext_corners.shape[1]-ii]
            # count += 1
            
        f.write('0\n')
                    
                 
                    
    # # Plot the geometry to make sure it looks as expected.
    # for ii in range(node_segments.shape[1]):
    #     if node_segments[0][ii] != 0 and node_segments[0][ii] != 0:
    #         plt.plot([all_the_nodes[0][int(node_segments[0][ii])-1], all_the_nodes[0][int(node_segments[1][ii])-1]], 
    #                  [all_the_nodes[1][int(node_segments[0][ii])-1], all_the_nodes[1][int(node_segments[1][ii])-1]],
    #                  color='r', linewidth=0.1)
    # plt.show()
        
                
                
    
    # Call Triangle to generate the mesh.
    # If triangle keeps bringing up errors, check the number of nodes and the
    # number of segments being reported in filename.poly
    # If we are running locally on a windows machine
    if sys.platform == 'win32':
        print("triangle -p -a{:.15f} -j -q30 -C {}".format(max_area, filename))
        t1 = time.time_ns()
        subprocess.run("triangle -p -a{:.15f} -j -q30 -C {}".format(max_area, filename))
        t2 = time.time_ns()
        print("Triangle completed in {:4.3g}s".format((t2-t1)*10**-9))
    # If we are running on BP (or any other linux machine). Commands to run are
    # incompatible between systems. Note also that if running on BP, the triangle
    # application needs to be compiled in the directory it will be run.
    elif sys.platform == 'linux':
        os.system("./triangle -p -a{:.15f} -j -q30 -C {}".format(max_area, filename))
        
        
        
    #%% Read all nodes and elements generated by Triangle
    
    # Get all node locations
    f = open('{}.1.node'.format(job_name_template), 'r')
    nodes = f.read().split('\n')
    node_lines = len(nodes)
    n_nodes = node_lines-3
    f.close()
    # Get all element locations
    f = open('{}.1.ele'.format(job_name_template), 'r')
    elements = f.read().split('\n')
    ele_lines = len(elements)
    n_elements = ele_lines-3
    f.close()
    
    # Store nodes in an easy to read format for determining if the node is in
    # the SRM. Disable if SRM is not in use, as this is very computationally
    # expensive and should be avoided if possible.
    node_block = np.zeros((2, n_nodes))
    for node in range(n_nodes):
        line = nodes[node+1].split('  ')
        line = [item for item in line if item != '']
        node_block[:, node] = [float(line[1]), float(line[2])]
        
    # Store elements in an easy to read format for SRM determination.
    element_block = np.zeros((3, n_elements))
    # Define array to store subset of elements which make up the SRM. Note that
    # an arbitrary large number has been used to initialise the array: after
    # all SRM elements have been obtained, the extras will need to be removed.
    SRM_elements = np.zeros((SRM_n_layers-1, 100000), dtype=int)
    # Subset of elements containing everything but the SRM elements. Again
    # initialised to be arbitrarily big.
    All_but_SRM_els = np.zeros(n_elements, dtype=int)
    # Ticker used to count how many elements are currently in each SRM sub-layer.
    SRM_els_ticker = np.zeros(SRM_n_layers-1, dtype=int)
    # Ticker for elements currently in use in All_but_SRM_els
    Not_ticker = 0
    # Ticker for current element.
    count = 0
    for element in range(n_elements):
        # Convert rows in elements to usable format.
        line = elements[element+1].split('  ')
        line = [item for item in line if item != '']
        # Store corner nodes. Element index is not stored as this is easily
        # obtained by how python indexes arrays (Triangle indexing starts from
        # 1, python indexing starts from 0).
        element_block[:, element] = [int(line[1]), int(line[2]), int(line[3])]
        
        # Check whether this element is in the SRM. Check whether all nodes are
        # contained by the SRM when doing this.
        SRM_ticker = 0
        el_pos = np.zeros(ext_corners[:, 0].shape)
        
        # If we found that this element was in a previously examined SRM,
        # skip the rest of the boundaries.
        if SRM_ticker != 0:
            break
        
        for node in range(3):
            # Nodes are recorded by Triangle from index 1, python refers from index 0.
            # Contents of element_block are recorded in Triangle notation.
            node_pos = node_block[:, int(element_block[node, element])-1]
        
            # Want this to be <= 0 if SRM is in -ve x-direction. Use machine
            # epsilon in case of equality.
            r = np.sqrt(node_pos[0]*node_pos[0] + node_pos[1]*node_pos[1])
            if r >= rad_at_SRM:
                SRM_ticker += 1
                el_pos[0] += node_pos[0] / 3
                el_pos[1] += node_pos[1] / 3
            # If any nodes are not, then this element is not in the SRM. Skip
            # the rest of the nodes in this element.
            else:
                break
                
        # If all nodes are in the SRM, check which sublayer the element is part of
        if SRM_ticker == 3:
            # The element must be at least in the 0th sublayer. As we may
            # have multiple SRMs, start indexing at -1 in case the element
            # is not in this SRM, but will be in another one we'll come
            # to on a later iteration.
            which_sublayer = -1
            # Check if it is in the next one.
            r = np.sqrt(el_pos[0]*el_pos[0] + el_pos[1]*el_pos[1])
            while r >= SRM_nodes[0][which_sublayer+1]:
                # If it is, update the current sublayer.
                which_sublayer += 1
            
            if which_sublayer > -1:
                # We now know which sublayer this element is in. Add it to that
                # sublayer. This element will not be in any more SRMs,
                # so we can stop checking.
                SRM_elements[which_sublayer, SRM_els_ticker[which_sublayer]] = int(line[0])
                SRM_els_ticker[which_sublayer] += 1
        elif SRM_ticker > 3:
            print("Node {} in too many SRMs - this message should never be printed to console.".format(int(line[0])))
        # If at least one node is not in the SRM, then the element is not in
        # the SRM. Store it in elsewhere.
        else:
            All_but_SRM_els[Not_ticker] = int(line[0])
            Not_ticker += 1
        
    # Reduce size of SRM_elements so that in 1st axis the array is as big as the
    # largest SRM subset. Note that there will still be zeros present in all
    # subsets but the largest. All zeros cannot be scrubbed because numpy
    # requires that arrays are rectangular, and scrubbing all would require rows
    # be different lengths.
    SRM_elements = SRM_elements[:, :np.max(SRM_els_ticker)]
    All_but_SRM_els = All_but_SRM_els[All_but_SRM_els != 0]
    
            
    
    #%% SRM sublayer parameters.
            
    # Get Young's modulus and Rayleigh Mass Damping values in SRM.
    SRM_modulus = bigE(0.01, SRM_thickness, lam_L, SRM_n_layers, modulus)
    SRM_CM = 2*np.pi*freq * bigX(SRM_n_layers, 3)
    
    
    
    #%% Start writing .inp files.
    
    print('{}'.format(job_name_template))
    
    # Work out approximate time taken for Abaqus to run a single job.
    Approx_time_per_el_per_step = 1.03e-7
    Est_runtime = n_elements * (step_time_length / time_step) * Approx_time_per_el_per_step
    if Est_runtime < 60:
        print('Est runtime per transmitter: {:.2f} s'.format(Est_runtime))
    elif Est_runtime < 3600:
        print('Est runtime per transmitter: {:.2f} m'.format(Est_runtime/60))
    else:
        print('Est runtime per transmitter: {:.2f} h'.format(Est_runtime/3600))   
    
    for el in range(num_els):
        
        # When corner is rounded, we're probably doing a diffraction study. As
        # such, we only need to activate one transmitter.
        if 'round' in settings['mesh']:
            if el+1 != 16:
                continue
        
        input_name = '{}_{}.inp'.format(job_name_template, el+1)
        
        # Open the input file and write the preamble
        with open(input_name, 'w') as i:
            i.write('*Heading\n')
            i.write('** Job name: {} Model name: Model-1\n'.format(job_name_template))
            i.write('*Preprint, echo=NO, model=NO, history=NO, contact=NO\n')
            i.write('**\n')
            i.write('** PART INSTANCE: Part-1-1\n')
            
            # Write all node locations
            i.write('*Node\n')
            for c1 in range(1,node_lines-2): #1 and -2 get rid of the extra bits
                line = list(filter(None, (nodes[c1]).split(' ')))
                i.write(', '.join(line[:3]) + '\n')
            
            #Write the element types
            i.write('*Element, type={}\n'.format(element_type))
            for c1 in range(1,ele_lines-2): #1 and -2 get rid of the extra bits
                line = filter(None, (elements[c1]).split(' '))
                i.write(', '.join(line) + '\n')    
                
            # Create node and element sets.
            
            i.write('*Nset, nset=All_nodes, generate\n')
            i.write('1, {}, 1\n'.format(n_nodes))
            
            i.write('*Elset, elset=All_elements, generate\n')
            i.write('1, {}, 1\n'.format(n_elements))
                
            # Create node sets corresponding to each element in the array
            for element in range(num_els):
                i.write('*Nset, nset=El{}\n'.format(element+1))
                for a in range(num_nodes_on_tx):
                    i.write('{},\n '.format(ext_corners.shape[1] + element*num_nodes_on_tx + a + 1))
                    
            upToMeas = ext_corners.shape[1]
            for layer in SRM_internal_nodes:
                upToMeas += layer.shape[1]
            i.write('*Nset, nset=MeasureSet\n')
            for a in range(meas_internal_nodes.shape[1]):
                i.write('{},\n '.format(upToMeas + a + 1))
                
            # Create element set for everything but the SRM, and each SRM
            # sublayer.
            i.write('*Elset, elset=All_but_SRM\n')
            for a in All_but_SRM_els:
                i.write('{},\n'.format(a))
            for sublayer in range(1, SRM_n_layers):
                i.write('*Elset, elset=SRM{}\n'.format(sublayer))
                for a in SRM_elements[sublayer-1]:
                    if a != 0:
                        i.write('{}, \n'.format(a))
            
            # Write the section definitions, matching regions of the mesh with
            # material definitions.
            i.write('*Solid Section, elset=All_but_SRM, material={}\n'.format(material_name))
            i.write(',\n')
            for sublayer in range(1, SRM_n_layers):
                i.write('*Solid Section, elset=SRM{}, material={}{}\n'.format(sublayer, material_name, sublayer))
                i.write(',\n')
                
            # Write the amplitude definition
            i.write('*System\n')
            i.write('*Amplitude, name=Amp-1\n')
            for a in amplitude:
                i.write('{}\n'.format(', '.join(map(str, a))))
            
            # Write the material definitions. Start with the bulk material.
            i.write('*Material, name={}\n'.format(material_name))
            i.write('*Density\n')
            i.write('{},\n'.format(density))
            i.write('*Elastic\n')
            i.write('{}, {}\n'.format(modulus, poisson))
            # Write the SRM materials. E and C_M will change as we get deeper
            # into the SRM.
            for sublayer in range(1, SRM_n_layers):
                i.write('*Material, name={}{}\n'.format(material_name, sublayer))
                i.write('*Density\n')
                i.write('{},\n'.format(density))
                i.write('*Elastic\n')
                i.write('{}, {}\n'.format(SRM_modulus[sublayer], poisson))
                i.write('*Damping, alpha={}\n'.format(SRM_CM[sublayer]))
            
            #Write the loading step
            i.write('*Step, name=Step-1, nlgeom=NO\n')
            i.write('*Dynamic, Explicit, direct user control\n')
            i.write('{}, {}\n'.format(time_step, step_time_length))
            
            # Write bulk viscosity, don't know if this is needed but it was
            # defined in the original copy of this script and it does no harm
            # so include it for now until it is determined to be needed or not.
            i.write('*Bulk Viscosity\n')
            i.write('0.0, 0.0\n')
            
            # Apply loads
            for element in range(num_els):
                i.write('*Cload, amplitude=Amp-1\n')
                i.write('El{}, 2, -10\n'.format(element+1)) # Node set, DoF, Magnitude
            
            # Outputs. Work through outputs requested in .yml file.
            if len(settings['output']) > 0:
                for outp in range(len(settings['output'])):
                    # Get output type and settings.
                    (thiskey, thisoutput) = list(settings['output'][outp].items())[0]
                    
                    if 'set' in thisoutput.keys():
                    
                        if thisoutput['type'] == 'field':
                            i.write('*Output, field, time interval={}\n'.format(thisoutput['t_int']))
                        elif thisoutput['type'] == 'history':
                            i.write('*Output, history, frequency={}\n'.format(thisoutput['freq']))
                        
                        if thisoutput['output'] == 'Element':
                            settype = 'elset'
                        elif thisoutput['output'] == 'Node':
                            settype = 'nset'
                    
                        if thisoutput['type'] == 'field':
                            i.write('*{} Output, {}=All_nodes\n'.format(thisoutput['output'], settype))
                        elif thisoutput['type'] == 'history':
                            i.write('*{} Output, {}=MeasureSet\n'.format(thisoutput['output'], settype))
                            
                        i.write('{},\n'.format(thiskey))
                        
                    elif 'sets' in thisoutput.keys():
                        elrange = [int(num) for num in thisoutput['sets'].split(', ')]
                        for element in range(elrange[0], elrange[1]+1):
                            if thisoutput['type'] == 'field':
                                i.write('*Output, field, time interval={}\n'.format(thisoutput['t_int']))
                            elif thisoutput['type'] == 'history':
                                i.write('*Output, history, frequency={}\n'.format(thisoutput['freq']))
                            
                            if thisoutput['output'] == 'Element':
                                settype = 'elset'
                            elif thisoutput['output'] == 'Node':
                                settype = 'nset'
                            
                            i.write('*{} Output, {}=El{}\n'.format(thisoutput['output'], settype, element))
                            i.write('{},\n'.format(thiskey))
            
            #End of step
            i.write('*End Step\n')
            
            print('Tx = {} written'.format(el+1))
    
    print('Input file written')
    
    
    
#%% Main run sequence

if __name__ == '__main__':
    # If python script is being run from the command line
    if len(sys.argv) > 1:
        yaml_name = '{}.yml'.format(sys.argv[1])
    # If script is being run from an editor
    else:
        yaml_name = 'O_1.yml'
        os.chdir('C:\\Users\\mc16535\\OneDrive - University of Bristol\\Documents\\Postgrad\\Coding\\Abaqus\\FMC Generation\\v8 - directivity')
    # Open and read .yml
    settings = read_settings(yaml_name)
    # Write inputs for all transmitters in FMC
    write_input(settings)