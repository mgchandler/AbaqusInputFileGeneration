function fn_orientation(fid, orient_name, a, b, rot_axis, rot_ang)
%USAGE
%   fn_orientation(fid, orient_name, a, b, rot_axis, rot_ang)
%SUMMARY
%   defines orientation
%INPUTS
%   fid - file ID
%   orient_name - name of orientation
%   a, b - coordinates of points (3 el vectors) to define orientation
%   (see abaqus manual). If left blank, these will be assumed to be
%   (1,0,0) and (0,1,0) which means orientation is same as global coords
%   rot_axis - local axis for coordinate rotation (1=x, 2=y, 3=z) default = 1
%   rot_ang - rotation angle in degrees about local axis, default = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = [fid, 1];
if isempty(a) & isempty(b)
    a = [1,0,0];
    b = [0,1,0];
end;
if isempty(rot_axis) & isempty(rot_ang)
    rot_axis = 1;
    rot_ang = 0;
end;
for ii = 1:length(fid)
    fprintf(fid(ii), ['*ORIENTATION, NAME=', orient_name, ', SYSTEM=RECTANGULAR\n']);
    fprintf(fid(ii), [repmat('%.5E, ', 1, 5), '%.5E\n'], a(:), b(:));
    fprintf(fid(ii), '%d, %.5E\n', round(rot_axis), rot_ang);
end;
return;
