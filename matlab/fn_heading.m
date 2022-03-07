function fn_heading(fid, heading_text)
%USAGE
%   fn_heading(fid, heading_text)
%SUMMARY
%   write header
%INPUTS
%   fid - file ID
%   heading_text - text to write in header, can include new line chars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = [fid, 1];
for ii = 1:length(fid)
    fprintf(fid(ii), '*HEADING\n');
    fprintf(fid(ii), heading_text);
end;
return;