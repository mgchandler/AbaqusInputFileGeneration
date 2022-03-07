function fn_elgen(fid, el_num1, ...
    num_els, nd_num_inc_in_row, el_num_inc_in_row, ...
    num_rows, nd_num_inc_by_row, el_num_inc_by_row, ...
    num_layers, nd_num_inc_by_layer, el_num_inc_by_layer, ...
    el_set, varargin)
%USAGE
%   fn_elgen(fid, el_num1, num_els, nd_num_inc_in_row, el_num_inc_in_row, num_rows, nd_num_inc_by_row, el_num_inc_by_row, num_layers, nd_num_inc_by_layer, el_num_inc_by_layer, el_set [, option_text])
%SUMMARY
%   element generation
%INPUTS
%   fid - file ID
%   num_els - number of elements in row
%   nd_num_inc_in_row - node number increment along row
%   el_num_inc_in_row - element number increment along row
%   num_rows - number of rows of elements
%   nd_num_inc_by_row - node number increment between rows
%   el_num_inc_by_row - element number increment between rows
%   num_layers - number of layers of elements
%   nd_num_inc_by_layer - node number increment between layers
%   el_num_inc_by_layer - element number increment between layers
%   el_set - name of element set (can be [])
%   [option_text] - optional text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 12
    option_text = varargin{1};
else
    option_text = '';
end;

fid = [fid, 1];
for ii = 1:length(fid)
    temp = '*ELGEN';
    if ~isempty(el_set)
        temp = [temp, ', ELSET=', el_set];
    end;
    if ~isempty(option_text)
        temp = [temp, ', ', option_text];
    end;
    temp = [temp, '\n'];
    fprintf(fid(ii), temp);
    %if the info fits on one line
    temp = [el_num1(:), ...
        num_els(:), nd_num_inc_in_row(:), el_num_inc_in_row(:), ...
        num_rows(:), nd_num_inc_by_row(:), el_num_inc_by_row(:), ...
        num_layers(:), nd_num_inc_by_layer(:), el_num_inc_by_layer(:)];
    fprintf(fid(ii), [repmat('%d, ', 1, size(temp, 2) - 1),'%d\n'], temp');
end;
return;
