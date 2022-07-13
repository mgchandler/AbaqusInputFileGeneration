function FMC_read(path, name, varargin)

disp(path)
disp(name)

cd(path)

% Get number of elements. Assume that filenames look like "txII-rxJJ.dat".
list = dir('tx*rx*.dat');
num_els_sq = length(list);
num_els = sqrt(num_els_sq);
assert(num_els == round(num_els), "FMC_read: incorrect number of files in directory.")

el = 0;

disp('Beginning data capture')

for tx = 1:num_els
    for rx = 1:num_els
        el = el+1;
        % Grab data
        fname = sprintf('tx%d-rx%d.dat', tx, rx);
        fID = fopen(fname);
        C = textscan(fID, '%f %f');
        t = C{1};
        d = C{2};
        fclose(fID);
        
        % Initialise arrays now that we know how big they are
        if and(tx == 1, rx == 1)
            arraylength = size(t, 1);
            time = zeros(arraylength, num_els^2);
            data = zeros(arraylength, num_els^2);
        end
        
        % Store data
        time(:, el) = t;
        data(:, el) = d;
    end
    fprintf('tx = %d complete\n', tx)
end

p = inputParser;
p.CaseSensitive = false;
addParameter(p, "freq", 5e6, @(x)validateattributes(x,{'numeric'},...
            {'nonempty','positive'}));
addParameter(p, "cycles", 5, @(x)validateattributes(x,{'numeric'},...
            {'nonempty','integer','nonnegative'}));
parse(p, varargin{:})

freq = p.Results.freq;
cycles = p.Results.cycles;


% Remove duplicate 0 term at start, and time shift
time = time(2:end, :) - .5*cycles/freq;
data = data(2:end, :);

% fft and ifft to get envelope and abs data.
fft_pts = 2^nextpow2(size(time, 1));

spec = fft(data, fft_pts, 1);
spec = spec(1:round(fft_pts/2), :);

sig = ifft(spec, fft_pts, 1);
sig = sig(1:size(time, 1), :);

data = sig;

% Leftover code for checking peak locations. May be used later.

% figure(1)
% imagesc(abs(sig))
% hold on
% locations = [];
% el = 0;
% for tx = 1:num_els
%     for rx = 1:num_els
%         el = el + 1;
%         
%         if tx == rx
%             [pks, locs] = findpeaks(squeeze(abs(sig(:, el))));
%             locations = [locations; locs];
%             scatter(ones(size(locs))*el, locs, 25, 'r', 'filled')
%         end
%     end
% end
% for ii = 1:num_els^2
%     [pks, locs] = findpeaks(squeeze(abs(sig(:, ii))));
%     locations = [locations; locs];
%     scatter(ones(size(locs))*ii, locs, 25, 'r', 'filled')
% end
% hold off

% locations(locations < 10150) = []; %5500, 10150
% locations(locations > 11000) = []; %6400, 11000

% figure(2)
% histogram(locations)

% disp('obs')
% disp(time(5975, 1))
% disp(time(10500, 1))
% disp('vis')
% disp(time(5875, 1))
% disp(time(10525, 1))

save(sprintf('%s.mat', name), 'time', 'data')

end