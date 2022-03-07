function fn_boundary_type(fid, node_set, bc_type)
%USAGE
%   fn_boundary_type(fid, node_set, bc_type)
%SUMMARY
%   set boundary conditions by type
%INPUTS
%   fid - file ID
%   node_set - name of node set
%   bc_type - type of boundary condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~iscell(node_set)
    node_set = cellstr(node_set);
    bc_type = cellstr(bc_type);
end;

fid = [fid, 1];
for ii = 1:length(fid)
    fprintf(fid(ii), '*BOUNDARY\n');
    for jj = 1:length(node_set)
        fprintf(fid(ii), [node_set{jj}, ', ', bc_type{jj}, '\n']);
    end;
end;
return;
