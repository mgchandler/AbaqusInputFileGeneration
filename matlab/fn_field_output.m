function fn_field_output(fid, node_set, number)
%USAGE
%   fn_field_output(fid, node_set, number)
%SUMMARY
%   write field output of nodal displacements to database
%INPUTS
%   fid - file ID
%   node_set - vector of node numbers
%   number - number of times field will be output will be generated over complete
%   marching step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = [fid, 1];
for ii = 1:length(fid)
    fprintf(fid(ii), '*OUTPUT, FIELD, NUMBER INTERVAL=%d\n', number);
    fprintf(fid(ii), ['*NODE OUTPUT, NSET=', node_set, '\n']);
    fprintf(fid(ii), 'U\n');
end;
return;