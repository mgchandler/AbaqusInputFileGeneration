function fn_start_time_march_step(fid, time_step, total_time)
%USAGE
%   fn_start_time_march_step(fid, time_step, total_time)
%SUMMARY
%   start time marching step
%INPUTS
%   fid - file ID
%   time_step - time step (if blank ABAQUS will calculate default value
%   total_time - total time to consider
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = [fid, 1];
for ii = 1:length(fid)
    fprintf(fid(ii), '*STEP, NLGEOM=NO\n');
    if isempty(time_step)
        fprintf(fid(ii), '*DYNAMIC, EXPLICIT, FIXED TIME INCREMENTATION\n');
        fprintf(fid(ii), ', %.5E\n', total_time);
    else
        fprintf(fid(ii), '*DYNAMIC, EXPLICIT, DIRECT USER CONTROL\n');
        fprintf(fid(ii), '%.5E, %.5E\n', time_step, total_time);
    end;
end;
return;