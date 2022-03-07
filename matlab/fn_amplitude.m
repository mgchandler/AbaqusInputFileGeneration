function fn_amplitude(fid, amp_name, time_step, centre_freq, number_of_cycles, window_type)
%USAGE
%   fn_amplitude(fid, amp_name, time_step, centre_freq, number_of_cycles,
%   window_type)
%SUMMARY
%   defines input signal amplitude curve (typically a toneburst or a gaussian
%   pulse
%INPUTS
%   fid - file ID
%   amp_name - name of amplitude curve
%   time_step - time step to use between points in curve
%   centre_freq - centre freq
%   number_of_cycles - number of cycles
%   window_type - 'hanning', 'gaussian', 'pulse'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
duration = number_of_cycles / centre_freq;
time = linspace(0, duration, ceil(duration / time_step));
centre_time = duration / 2;
sig = sin(2 * pi * centre_freq * (time - centre_time));
switch window_type
    case 'hanning'
        win = 0.5 * (1 + cos(2 * pi * centre_freq * (time - centre_time) / number_of_cycles));
        win(find(abs(time - centre_time) > duration / 2)) = 0;
        sig = sig .* win;
    case 'gaussian'
        error('Not implemented yet');
        sig = sig .* win;
    case 'pulse'
        error('Not implemented yet');
end;
data = [time(:), sig(:)];
fid = [fid, 1];
for ii = 1:length(fid)
    fprintf(fid(ii), ['*AMPLITUDE, NAME=', amp_name, '\n']);
    for jj = 1:floor(size(data, 1) / 4)
        temp = data((jj-1) * 4 + 1 : jj * 4, :);
        fprintf(fid(ii), [repmat('%.5E, ', 1, prod(size(temp)) - 1), '%.5E\n'], temp');
        if fid(ii) == 1 & jj > 2
            fprintf(fid(ii), '... some amplitude lines suppressed on screen\n');
            break;
        end;
    end;
    if rem(size(data, 1), 4)
        temp = data(end - rem(size(data, 1), 4): end, :);
        fprintf(fid(ii), [repmat('%.5E, ', 1, prod(size(temp)) - 1), '%.5E\n'], temp');
    end;
end;
return;
