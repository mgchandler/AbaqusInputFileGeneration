function fn_create_2d_rect_region(fid, xmin, xmax, ymin, ymax, el_size, name)
%USAGE
%   fn_create_2d_rect_region(fid, xmin, xmax, ymin, ymax, el_size, name)
%SUMMARY
%   creates a 2D rectangular region of plane strain elements with following
%   node sets:
%       NSET_name_ALL
%       NSET_name_EDGE_XMIN
%       NSET_name_EDGE_XMAX
%       NSET_name_EDGE_YMIN
%       NSET_name_EDGE_YMAX
%   and element sets:
%       ELSET_name_ALL
%INPUTS
%   fid - file ID
%   xmin, xmax, ymin, ymax - positions of edges
%   el_size - maximum element size (elements will be less than
%   or equal to this in size in both dimensions)
%   name - name used in node and element set names (can be [])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base_nset_name = ['NSET_', name];
base_elset_name = ['ELSET_', name];

x_els = ceil((xmax - xmin) / el_size);
y_els = ceil((ymax - ymin) / el_size);
x_el_size = (xmax - xmin) / x_els;
y_el_size = (ymax - ymin) / y_els;
x_nodes = x_els + 1;
y_nodes = y_els + 1;

%bottom row
index1 = 1;
index2 = index1 + x_nodes - 1;
node_num = [index1, index2]';
node_pos = [[xmin, ymin]; [xmax, ymin]];
node_set = [];
fn_node(fid, node_num, node_pos, node_set)
node_set = [base_nset_name, '_EDGE_YMIN'];
increment = 1;
fn_ngen(fid, index1, index2, increment, node_set);

%copy row to generate top row
old_node_set = [base_nset_name, '_EDGE_YMIN'];
option_text = ['NEW SET=', base_nset_name, '_EDGE_YMAX'];
change_number = x_nodes * (y_nodes - 1);
dx = 0.0;
dy = ymax - ymin;
dz = 0.0;
fn_ncopy_shift(fid, change_number, old_node_set, dx, dy, dz, option_text);

%fill in gap
node_set1 = [base_nset_name, '_EDGE_YMIN'];
node_set2 = [base_nset_name, '_EDGE_YMAX'];
node_set = [base_nset_name, '_ALL'];
intervals_along_line = y_nodes - 1;
increment_between_lines = x_nodes;
fn_nfill(fid, node_set1, node_set2, intervals_along_line, increment_between_lines, node_set);

%name other node set
node_set = [base_nset_name, '_EDGE_XMIN'];
node_num1 = 1;
node_num2 = x_nodes * (y_nodes - 1) + 1;
increment = x_nodes;
fn_nset(fid, node_set, node_num1, node_num2, increment);

node_set = [base_nset_name, '_EDGE_XMAX'];
node_num1 = x_nodes;
node_num2 = x_nodes * y_nodes;
increment = x_nodes;
fn_nset(fid, node_set, node_num1, node_num2, increment);

%do the first element
element_type = 'CPE4R';
element_num = 1;
node_nums = [1, 2, x_nodes + 2, x_nodes + 1];
fn_element(fid, element_type, element_num, node_nums, [])

%fill in rest of mesh
el_num1 = 1;
num_els = x_els;
nd_num_inc_in_row = 1;
el_num_inc_in_row = 1;
num_rows = y_els;
nd_num_inc_by_row = x_nodes;
el_num_inc_by_row = x_els;
num_layers = 1;
nd_num_inc_by_layer = 1;
el_num_inc_by_layer = 1;
el_set = [base_elset_name, '_ALL'];
fn_elgen(fid, el_num1, num_els, nd_num_inc_in_row, el_num_inc_in_row, num_rows, nd_num_inc_by_row, el_num_inc_by_row, num_layers, nd_num_inc_by_layer, el_num_inc_by_layer, el_set);
return;