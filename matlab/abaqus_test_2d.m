clc;
close all;
clear;

%test program to write abaqus input file
run_job = 0; 
input_fname = 'test_2d.inp';

x_size = 0.1;
y_size = 0.001;
el_size = 0.001/10;
% x_size = 0.004;
% y_size = 0.003;
% el_size = 0.001;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_els = ceil(x_size / el_size);
y_els = ceil(y_size / el_size);
x_nodes = x_els + 1;
y_nodes = y_els + 1;

fid = fopen(input_fname, 'wt');
if fid < 0
    error('Could not create input file - previous job may not have terminated');
    return;
end;
heading_text = sprintf('This is a test.\nPlate %.5f x %.5f mm\n%.5f mm elements\n', x_size * 1e3, y_size * 1e3, el_size * 1e3);
fn_heading(fid, heading_text);

fn_create_2d_rect_region(fid, 0, x_size, 0, y_size, el_size, 'PDW');

%define material
matl_name = 'STEEL';
matl_props.density = 7800;
matl_props.youngs_modulus = 207e9;
matl_props.poissons_ratio = 0.33;
%cross ply
% matl_name = 'CROSS_PLY_COMPOSITE';
% matl_props.density = 1560;
% matl_props.stiffness_matrix = [ ...
%     [64.24, 5.6, 7.73, 0, 0, 0]; ...
%     [5.6, 70.87, 8.39, 0, 0, 0]; ...
%     [7.73, 8.39, 12.15, 0, 0, 0]; ...
%     [0, 0, 0, 2.97, 0, 0]; ...
%     [0, 0, 0, 0, 3.06, 0]; ...
%     [0, 0, 0, 0, 0, 4.7]] * 1e9;
fn_define_material(fid, matl_name, matl_props);

%define section
el_set = 'ELSET_PDW_ALL';
orient_name = 'ORIENT_GLOBAL';
fn_define_solid_section(fid, el_set, matl_name, orient_name);

%define orientation
fn_orientation(fid, orient_name, [], [], [], []);

%define amplitude of input
amp_name = '5_CYCLE_1_MHZ_HANNING';
time_step = 1 / 1e6 / 10;
centre_freq = 1e6;
number_of_cycles = 5;
window_type = 'hanning';
fn_amplitude(fid, amp_name, time_step, centre_freq, number_of_cycles, window_type);

%time-marching step
time_step = [];
total_time = 10e-6;
fn_start_time_march_step(fid, time_step, total_time);

%boundary conditions
node_set = {'NSET_PDW_EDGE_XMIN', 'NSET_PDW_EDGE_XMAX'};
bc_type = {'XSYMM', 'XSYMM'};
fn_boundary_type(fid, node_set, bc_type);

%apply load
amp_name = '5_CYCLE_1_MHZ_HANNING';
node_num_or_set = x_nodes * (y_nodes - 1) + (x_nodes + 1) / 2;
dof = 2;
amp = 1;
fn_cload(fid, amp_name, node_num_or_set, dof, amp);

%specify output
node_set = 'NSET_PDW_ALL';
number = 10;
fn_field_output(fid, node_set, number);

%end step
fn_end_step(fid);

fclose(fid);

if run_job == 1
    job_fname = fn_run_abaqus(input_fname)
end;

if run_job == 2
    job_fname = fn_check_abaqus(input_fname)
end;