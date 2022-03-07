function job_fname = fn_check_abaqus(input_fname, varargin)
if nargin < 2
    [pathstr, job_fname, ext, version] = fileparts(input_fname);
else
    job_fname = varargin{1};
end;
str = ['!c:\abaqus\commands\abaqus', ...
    ' job=',job_fname, ...
    ' input=',input_fname, ...
    ' datacheck'];
eval(str);
return;
