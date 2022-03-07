function fn_elset(fid, el_set, el_num1, el_num2, varargin)
%USAGE
%   fn_elset(fid, el_set, el_num1, el_num2 [, increment, option_text])
%SUMMARY
%   create set of elements from list (ABAQUS elset GENERATE option)
%INPUTS
%   fid - file ID
%   el_set - name of element set
%   node_num1 - number(s) of first node in group(s)
%   node_num2 - number(s) of last node in group(s)
%   [increment = 1] - increment between node numbers in set
%   [option_text] - optional text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 4
    increment = varargin{1};
else
    increment = zeros(length(el_num1),1);
end;
if nargin > 5
    option_text = varargin{2};
else
    option_text = '';
end;

fid = [fid, 1];
for ii = 1:length(fid)
    if isempty(option_text)
        fprintf(fid(ii), ['*ELSET, ELSET=',el_set,', GENERATE\n']);
    else
        fprintf(fid(ii), ['*ELSET, ELSET=',el_set,', GENERATE, ', option_text, '\n']);
    end;
    for jj = 1:length(el_num1)
        if ~increment(jj)
            fprintf(fid(ii), ['%d, %d\n'], el_num1, el_num2);
        else
            fprintf(fid(ii), ['%d, %d, %d\n'], el_num1, el_num2, increment(jj));
        end;
    end;
end;
return;
