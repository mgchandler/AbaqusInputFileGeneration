function fn_create_2d_rect_region(fid, xmin, xmax, ymin, ymax, zmin, zmax, el_size, name)
%USAGE
%   fn_create_3d_rect_region(fid, xmin, xmax, ymin, ymax, zmin, zmax, el_size, name)
%SUMMARY
%   creates a 3D rectangular region of brick elements with following
%   node sets:
%       NSET_name_ALL
%       NSET_name_FACE_XMIN
%       NSET_name_FACE_XMAX
%       NSET_name_FACE_YMIN
%       NSET_name_FACE_YMAX
%       NSET_name_FACE_ZMIN
%       NSET_name_FACE_ZMAX
%       NSET_name_EDGE_XMIN_YMIN
%       NSET_name_EDGE_XMIN_YMAX
%       NSET_name_EDGE_XMAX_YMIN
%       NSET_name_EDGE_XMAX_YMAX
%       NSET_name_EDGE_XMIN_ZMIN
%       NSET_name_EDGE_XMIN_ZMAX
%       NSET_name_EDGE_XMAX_ZMIN
%       NSET_name_EDGE_XMAX_ZMAX
%       NSET_name_EDGE_YMIN_ZMIN
%       NSET_name_EDGE_YMIN_ZMAX
%       NSET_name_EDGE_YMAX_ZMIN
%       NSET_name_EDGE_YMAX_ZMAX
%   and element sets:
%       ELSET_name_ALL
%INPUTS
%   fid - file ID
%   xmin, xmax, ymin, ymax, zmin, max - positions of edges
%   el_size - maximum element size (elements will be less than
%   or equal to this in size in both dimensions)
%   name - name used in node and element set names (can be [])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base_nset_name = ['NSET_', name];
base_elset_name = ['ELSET_', name];

%basic size of model
x_els = ceil((xmax - xmin) / el_size);
y_els = ceil((ymax - ymin) / el_size);
z_els = ceil((zmax - zmin) / el_size);
x_nodes = x_els + 1;
y_nodes = y_els + 1;
z_nodes = z_els + 1;

%various node indices to define corners of model
dx = 1;
dy = x_nodes;
dz = x_nodes * y_nodes;
xminyminzmin = 1;
xmaxyminzmin = x_nodes;
xminymaxzmin = x_nodes * (y_nodes - 1) + 1;
xmaxymaxzmin = x_nodes * y_nodes;
xminyminzmax = 1 + x_nodes * y_nodes * (z_nodes - 1);
xmaxyminzmax = x_nodes + x_nodes * y_nodes * (z_nodes - 1);
xminymaxzmax = x_nodes * (y_nodes - 1) + 1 + x_nodes * y_nodes * (z_nodes - 1);
xmaxymaxzmax = x_nodes * y_nodes + x_nodes * y_nodes * (z_nodes - 1);

%generate ymin, zmin row
fn_node(fid, [xminyminzmin, xmaxyminzmin]', [[xmin, ymin, zmin]; [xmax, ymin, zmin]], [])
fn_ngen(fid, xminyminzmin, xmaxyminzmin, dx, [base_nset_name, '_EDGE_YMIN_ZMIN']);

%copy ymin, zmin row to generate ymax, zmin row
fn_ncopy_shift(fid, x_nodes * (y_nodes - 1), [base_nset_name, '_EDGE_YMIN_ZMIN'], 0.0, ymax - ymin, 0.0, ['NEW SET=', base_nset_name, '_EDGE_YMAX_ZMIN']);

%fill in gap to create zmin surface
fn_nfill(fid, [base_nset_name, '_EDGE_YMIN_ZMIN'], [base_nset_name, '_EDGE_YMAX_ZMIN'], y_nodes - 1, x_nodes, [base_nset_name, '_FACE_ZMIN']);

%copy zmin face to create zmax face
fn_ncopy_shift(fid, x_nodes * y_nodes * (z_nodes - 1), [base_nset_name, '_FACE_ZMIN'], 0.0, 0.0, zmax - zmin, ['NEW SET=', base_nset_name, '_FACE_ZMAX']);

%fill in gap to complete node definition
node_set1 = [base_nset_name, '_FACE_ZMIN'];
node_set2 = [base_nset_name, '_FACE_ZMAX'];
node_set = [base_nset_name, '_ALL'];
intervals_along_line = z_nodes - 1;
increment_between_lines = x_nodes * y_nodes;
fn_nfill(fid, node_set1, node_set2, intervals_along_line, increment_between_lines, node_set);

%do the first element
element_type = 'C3D8R';
element_num = 1;
node_nums = [1, 2, x_nodes + 2, x_nodes + 1];
node_nums = [node_nums, node_nums + x_nodes * y_nodes];
fn_element(fid, element_type, element_num, node_nums, [])

%fill in rest of mesh
el_num1 = 1;
num_els = x_els;
nd_num_inc_in_row = 1;
el_num_inc_in_row = 1;
num_rows = y_els;
nd_num_inc_by_row = x_nodes;
el_num_inc_by_row = x_els;
num_layers = z_els;
nd_num_inc_by_layer = x_nodes * y_nodes;
el_num_inc_by_layer = x_els * y_els;
el_set = [base_elset_name, '_ALL'];
fn_elgen(fid, el_num1, num_els, nd_num_inc_in_row, el_num_inc_in_row, num_rows, nd_num_inc_by_row, el_num_inc_by_row, num_layers, nd_num_inc_by_layer, el_num_inc_by_layer, el_set);

%create sets for other edges
fn_nset(fid, [base_nset_name, '_EDGE_XMIN_YMIN'], xminyminzmin, xminyminzmax, dz);
fn_nset(fid, [base_nset_name, '_EDGE_XMAX_YMIN'], xmaxyminzmin, xmaxyminzmax, dz);

fn_nset(fid, [base_nset_name, '_EDGE_XMIN_YMAX'], xminymaxzmin, xminymaxzmax, dz);
fn_nset(fid, [base_nset_name, '_EDGE_XMAX_YMAX'], xmaxymaxzmin, xmaxymaxzmax, dz);

fn_nset(fid, [base_nset_name, '_EDGE_XMIN_ZMIN'], xminyminzmin, xminymaxzmin, dy);
fn_nset(fid, [base_nset_name, '_EDGE_XMAX_ZMIN'], xmaxyminzmin, xmaxymaxzmin, dy);

fn_nset(fid, [base_nset_name, '_EDGE_XMIN_ZMAX'], xminyminzmax, xminymaxzmax, dy);
fn_nset(fid, [base_nset_name, '_EDGE_XMAX_ZMAX'], xmaxyminzmax, xmaxymaxzmax, dy);

fn_nset(fid, [base_nset_name, '_EDGE_YMIN_ZMAX'], xminyminzmax, xmaxyminzmax, dx);
fn_nset(fid, [base_nset_name, '_EDGE_YMAX_ZMAX'], xminymaxzmax, xmaxymaxzmax, dx);

%create sets for the other faces
fn_nset(fid, [base_nset_name, '_FACE_XMIN'], [xminyminzmin:dy:xminymaxzmin]', [xminyminzmax:dy:xminymaxzmax]', dz * ones(y_nodes, 1));
fn_nset(fid, [base_nset_name, '_FACE_XMAX'], [xmaxyminzmin:dy:xmaxymaxzmin]', [xmaxyminzmax:dy:xmaxymaxzmax]', dz * ones(y_nodes, 1));

fn_nset(fid, [base_nset_name, '_FACE_YMIN'], [xminyminzmin:dx:xmaxyminzmin]', [xminyminzmax:dx:xmaxyminzmax]', dz * ones(x_nodes, 1));
fn_nset(fid, [base_nset_name, '_FACE_YMAX'], [xminymaxzmin:dx:xmaxymaxzmin]', [xminymaxzmax:dx:xmaxymaxzmax]', dz * ones(x_nodes, 1));
return;