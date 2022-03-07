function fn_node(fid, node_num, node_pos, node_set, varargin)
%USAGE
%   fn_node(fid, node_num, node_pos, node_set [, option_text])
%SUMMARY
%   write node details
%INPUTS
%   fid - file ID
%   node_num - vector of node numbers
%   node_pos - vector of node positions. must have number of rows equal to
%   length(node_num) and columns for x, y and z positions as required
%   node_set - name of node set (can be [])
%   [option_text] - optional text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 4
    option_text = varargin{1};
else
    option_text = '';
end;
fid = [fid, 1];
for ii = 1:length(fid)
    temp = '*NODE';
    if ~isempty(node_set)
        temp = [temp, ', NSET=', node_set];
    end;
    if ~isempty(option_text)
        temp = [temp, ', ', option_text];
    end;
    temp = [temp, '\n'];
    fprintf(fid(ii), temp);
    temp = [node_num(:), node_pos];
    fprintf(fid(ii), ['%d, ',repmat('%.5E, ', 1, size(node_pos, 2) - 1),'%.5E\n'], temp');
end;
return;