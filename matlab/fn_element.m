function fn_element(fid, element_type, element_num, node_nums, el_set, varargin)
%USAGE
%   fn_element(fid, element_type, element_num, node_nums, varargin)
%SUMMARY
%   create elements at specified nodes
%INPUTS
%   fid - file ID
%   element_type - name of element type (e.g. 'CPE4R');
%   element_num - vector of element numbers
%   node_nums - matrix of node numbers (must have number of rows =
%   length(element_num)
%   el_set - name of element set (can be [])
%   [option_text] - optional text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 5
    option_text = varargin{1};
else
    option_text = '';
end;

fid = [fid, 1];
for ii = 1:length(fid)
    temp = ['*ELEMENT, TYPE=', element_type];
    if ~isempty(el_set)
        temp = [temp, ', ELSET=', el_set];
    end;
    if ~isempty(option_text)
        temp = [temp, ', ', option_text];
    end;
    temp = [temp, '\n'];
    fprintf(fid(ii), temp);
    if size(node_nums < 15)
        %if the info fits on one line
        temp = [element_num(:), node_nums];
        fprintf(fid(ii), ['%d, ',repmat('%d, ', 1, size(node_nums, 2) - 1),'%d\n'], temp');
    else
        %if it doesn't ...
        error('Not implemented yet')
    end;
end;
return;
