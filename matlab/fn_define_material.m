function fn_define_material(fid, matl_name, matl_props)
%USAGE
%   fn_define_material(fid, matl_name, matl_props)
%SUMMARY
%   defines material
%INPUTS
%   fid - file ID
%   matl_name - name of material
%   matl_props - structured variable. Must contain field matl_props.density
%   and may contain EITHER
%       matl_props.youngs_modulus AND matl_props.poissons_ratio
%OR     matl_props.stiffness_matrix - 6x6 stiffness matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = [fid, 1];
for ii = 1:length(fid)
    fprintf(fid(ii), ['*MATERIAL, NAME=', matl_name, '\n']);
    if isfield(matl_props, 'youngs_modulus')
        fprintf(fid(ii), '*ELASTIC, TYPE=ISOTROPIC\n');
        fprintf(fid(ii), '%.5E, %.5E\n', matl_props.youngs_modulus, matl_props.poissons_ratio);
    end;
    if isfield(matl_props, 'stiffness_matrix')
        fprintf(fid(ii), '*ELASTIC, TYPE=ANISOTROPIC\n');
        c = matl_props.stiffness_matrix;
        temp = [c(1, 1); c(1:2, 2); c(1:3, 3); c(1:4, 4); c(1:5, 5); c(1:6, 6)];
        fprintf(fid(ii), [repmat('%.5E, ', 1, 7), '%.5E\n'], temp(1:8));
        fprintf(fid(ii), [repmat('%.5E, ', 1, 7), '%.5E\n'], temp(9:16));
        fprintf(fid(ii), [repmat('%.5E, ', 1, 4), '%.5E\n'], temp(17:21));
%         error('Not implemented yet!');
    end;
    fprintf(fid(ii), '*DENSITY\n');
    fprintf(fid(ii), '%.5E\n', matl_props.density);
end;
return;
