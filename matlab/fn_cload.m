function fn_cload(fid, amp_name, node_num_or_set, dof, amp);
%USAGE
%   fn_cload(fid, amp_name, node_num, dof, amp)
%SUMMARY
%   apply concentrated load at specified node(s)
%INPUTS
%   fid - file ID
%   amp_name - name of amplitude data
%   node_num_or_set - vector of node numbers or cell array of node set
%   names
%   dof - vector of DOFs in which load is applied
%   amp - vector of amps of load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = [fid, 1];
if ~isnumeric(node_num_or_set) & length(node_num_or_set) == 1
    node_num_or_set = cellstr(node_num_or_set);
end;
for ii = 1:length(fid)
    fprintf(fid(ii), ['*CLOAD, AMPLITUDE=', amp_name, '\n']);
    for jj = 1:length(node_num_or_set)
        if iscell(node_num_or_set)
            temp = node_num_or_set{jj};
        else
            temp = sprintf('%d', node_num_or_set(jj));
        end;
        fprintf(fid(ii), [temp, ', %d, %.5E\n'], dof(jj), amp(jj));
    end;
end;
return;