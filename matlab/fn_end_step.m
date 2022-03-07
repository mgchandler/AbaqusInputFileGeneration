function fn_end_step(fid)
%USAGE
%   fn_end_step(fid)
%SUMMARY
%   ends step
%INPUTS
%   fid - file ID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = [fid, 1];
for ii = 1:length(fid)
    fprintf(fid(ii), '*END STEP\n');
end;
return;
