function fn_define_solid_section(fid, el_set, matl_name, orient_name)
%USAGE
%   fn_define_solid_section(fid, el_set, matl)
%SUMMARY
%   define solid section using normal controls for explicit
%INPUTS
%   fid - file ID
%   el_set - name of element set that makes up section
%   matl_name - material name
%   orient_name - orientation name (can be [])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = [fid, 1];
for ii = 1:length(fid)
    fprintf(fid(ii), '*SECTION CONTROLS, NAME=HG, HOURGLASS=STIFFNESS\n');
    fprintf(fid(ii), ['*SOLID SECTION, ELSET=', el_set,', MATERIAL=', matl_name,', CONTROLS=HG']);
    if isempty(orient_name)
        fprintf(fid(ii), '\n');
    else
        fprintf(fid(ii), [', ORIENTATION=', orient_name, '\n']);
    end;
end;
return;
