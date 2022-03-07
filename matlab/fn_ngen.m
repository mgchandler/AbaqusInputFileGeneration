function fn_ngen(fid, index1, index2, increment, node_set, varargin)
%USAGE
%   fn_ngen(fid, index1, index2, node_set, increment [, option_text])
%SUMMARY
%   write node details
%INPUTS
%   fid - file ID
%   index1 - index of first node(s)
%   index2 - index of first node(s)
%   node_set - name of node set (can be [])
%   increment - increment between nodes
%   [option_text] - optional text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 5
    option_text = varargin{1};
else
    option_text = '';
end;
fid = [fid, 1];
for ii = 1:length(fid)
    temp = '*NGEN';
    if ~isempty(node_set)
        temp = [temp, ', NSET=', node_set];
    end;
    if ~isempty(option_text)
        temp = [temp, ', ', option_text];
    end;
    temp = [temp, '\n'];
    fprintf(fid(ii), temp);
    if ~increment
        temp = [index1(:), index2(:)];
        fprintf(fid(ii), ['%d, %d\n'], temp');
    else
        temp = [index1(:), index2(:), increment(:)];
        fprintf(fid(ii), ['%d, %d, %d\n'], temp');
    end;
end;
return;