function FMC_read(path, name)

disp(path)
disp(name)

cd(path)

num_els = 32;

el = 0;

disp('Beginning data capture')

for tx = 1:num_els
    for rx = 1:num_els
        el = el+1;
        % Grab data
        fname = sprintf('tx%d-rx%d.dat', tx, rx);
        [t, d] = textread(fname, '%f %f');
        
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

% Remove duplicate 0 term at start, and time shift
time = time(2:end, :) - 5.0e-7;
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