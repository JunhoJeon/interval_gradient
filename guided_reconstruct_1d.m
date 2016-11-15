function [r, a, b, p] = guided_reconstruct_1d(g, gx, wx, sigma_s, epsi)

ss = sigma_s;

% reconstruct p from gx
p = g;
cgx = cumsum(gx, 2);
p(:, 2:end, :) = bsxfun(@plus, p(:, 1, :), cgx(:, 1:end-1, :));

q = g;

% 1d guided filter
mean_w = wx.^2; %boxfilter1d(wx.^2, ss);
mean_p = boxfilter1d(p, ss);
mean_pq = boxfilter1d(p.*q, ss);
mean_q = boxfilter1d(q, ss);
% mean_p = gaussian1d(p, ss);
% mean_pq = gaussian1d(p.*q, ss);
% mean_q = gaussian1d(q, ss);
cov_pq = max(0, mean_pq - mean_p.*mean_q);
mean_p2 = boxfilter1d(p.*p, ss);
% mean_p2 = gaussian1d(p.*p, ss);
var_p = max(0, mean_p2 - mean_p.^2);

epsi = epsi .* mean_w;

a = cov_pq ./ (var_p + epsi);
max_a = min(1, max(a, [], 3));
max_a = repmat(max_a, [1, 1, size(g, 3)]);
a = max(max_a, a);
% a(idx) = min(1, max_a(idx));
% a = repmat(a, [1 1 size(p, 3)]);

b = mean_q - a.*mean_p;

mean_a = boxfilter1d(a, ss);
mean_b = boxfilter1d(b, ss);
% mean_a = gaussian1d(a, ss);
% mean_b = gaussian1d(b, ss);

r = mean_a.*p + mean_b;

end
