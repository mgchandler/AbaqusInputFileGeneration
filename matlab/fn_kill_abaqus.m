function fn_kill_abaqus(job_fname)
    str = ['!c:\abaqus\commands\abaqus terminate', ...
        ' job=',job_fname];
    eval(str);
return;
