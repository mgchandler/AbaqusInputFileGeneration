function fn_nfill(fid, node_set1, node_set2, intervals_along_line, increment_between_lines, node_set, varargin)
%USAGE
%   fn_nfill(fid, node_set1, node_set2, intervals_along_line, increment_between_lines [, option_text])
%SUMMARY
%   fills in grid of nodes between two lines (defined by node_set names)
%INPUTS
%   fid - file ID
%   node_set1 - name of first node set(s) (cell array of strings)
%   node_set2 - name of second node set(s) (cell array of strings)
%   intervals_along_line - number of intervals along each generated line (must
%   be same length as node_set)
%   increment_between_lines - increment in node number between lines
%   node_set - name of generated set (can be []);
%   [option_text] - optional text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 6
    option_text = varargin{1};
else
    option_text = '';
end;

if ~iscell(node_set1)
    node_set1 = cellstr(node_set1);
end;
if ~iscell(node_set2)
    node_set2 = cellstr(node_set2);
end;
fid = [fid, 1];
for ii = 1:length(fid)
    temp = '*NFILL';
    if ~isempty(node_set)
        temp = [temp, ', NSET=', node_set];
    end;
    if ~isempty(option_text)
        temp = [temp, ', ', option_text];
    end;
    temp = [temp, '\n'];
    fprintf(fid(ii), temp);
    for jj = 1:length(node_set1)
        if ~increment_between_lines(jj)
            fprintf(fid(ii), [node_set1{jj}, ', ', node_set2{jj}, ', %d\n'], intervals_along_line(jj));
        else
            fprintf(fid(ii), [node_set1{jj}, ', ', node_set2{jj}, ', %d, %d\n'], intervals_along_line(jj), increment_between_lines(jj));
        end;
    end;
end;
return;