dx = 0.0012634024449781562 / 30;
topwallhalflen = 5e-3;
pitch = 1e-3;
decay_len = 0.0e-3;
N = round(2*topwallhalflen/dx);

xs = linspace(-topwallhalflen, topwallhalflen, N)';
xs = xs(2:end-1);
N = N-2;

field = zeros(N,4);
field(:, 1) = xs;

half_element_idxs = round(pitch/dx/2);

half_N = round(N/2);
n=round(decay_len / dx);
cosine = (1 + cos(pi * [0:n-1]/n)') / 2;

field(half_N - half_element_idxs : half_N + half_element_idxs, 4) = 1;

field(half_N - half_element_idxs - (n-1): half_N - half_element_idxs, 4) = flip(cosine);
field(half_N + half_element_idxs : half_N + half_element_idxs + (n-1), 4) = cosine;

