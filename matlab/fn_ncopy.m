function fn_ncopy_shift(fid, change_number, old_node_set, new_node_set, dx, dy, dz, varargin)
%USAGE
%   fn_ncopy_shift(fid, node_num, node_pos, node_set [, option_text])
%SUMMARY
%   copy node details by translation
%INPUTS
%   fid - file ID
%   change_number - increment in node number between copies
%   old_node_set - name of original node set
%   new_node_set - name of new node set
%   dx, dy, dz - translation to apply
%   [option_text] - optional text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 7
    option_text = varargin{1};
else
    option_text = '';
end;
fid = [fid, 1];
for ii = 1:length(fid)
    temp = sprintf(['*NCOPY, OLD SET=', old_node_set, ', CHANGE NUMBER=%d, SHIFT'], change_number);
    if ~isempty(node_set)
        temp = [temp, ', NEW SET=', node_set];
    end;
    if ~isempty(option_text)
        temp = [temp, ', ', option_text];
    end;
    temp = [temp, '\n'];
    fprintf(fid(ii), temp);
    fprintf(fid, '%.5E, %.5E, %.5E\n', dx, dy, dz);
    
end;
return;