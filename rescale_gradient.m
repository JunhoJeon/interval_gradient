function [r_gx, wx] = rescale_gradient(I, ss, grad_type, do_normalize, I0)

if ~exist('do_normalize', 'var'), do_normalize = false; end

% ss = ss*sqrt(2);
fr = ceil(3*ss);

% gaussian patch graidient
switch grad_type
  case 'gaussian'
    kl = -fr+1:0;
    kl = exp(-0.5*(kl/ss).^2);
    kl = kl / sum(kl);
    k = [0 -kl kl(end:-1:1)];
    
    kg = zeros(1, length(kl)*2-1);
    for i = 0:length(kl)-1
      kg(1+i:end-i) = kg(1+i:end-i) + kl(i+1);
    end
%     kg = kg / sum(kg) * 2;
  case 'box'
    frb = ceil(2*ss);
    kl = ones(1, frb);
    k = [0 -kl kl];
    k = k / sum(abs(k)) * 2;
  case 'dog'
    x = -fr:fr;
    k = exp(-0.5*(x/ss).^2);
    k = x .* k / (ss.^2);
    k = k / sum(abs(k)) * 2;
    k(2:fr+1) = k(1:fr);
    k(1) = 0;
    
    kg = zeros(1, length(k)-2);
    for i = 0:floor(length(kg)/2)
      kg(1+i:end-i) = kg(1+i:end-i) - k(2+i);
    end
%     kg = kg / sum(kg);
  case 'rtv'
    [r_gx, wx] = rescale_gradient_rtv(I, ss, do_normalize);
    return;
  otherwise
    disp('unknown interval gradient type, using gaussian.');
    kl = -fr+1:0;
    kl = exp(-0.5*(kl/ss).^2);
    kl = kl / sum(kl);
    k = [0 -kl kl(end:-1:1)];
end

ky = fspecial('gaussian', [2*fr+1 1], ss);

gx = imfilter(I, [0 -1 1], 'replicate');
px = imfilter(I, k, 'replicate');
% px = imfilter(px, ky, 'symmetric');

% dx = imfilter(abs(gx), kg, 'symmetric');
% px = imfilter(gx, kg, 'symmetric');

nch = size(I, 3);

[~, midx] = max(abs(gx), [], 3);

[hh, ww, ~] = size(gx);
pidx = reshape(1:hh*ww, [hh ww]);
midx = (midx-1)*hh*ww + pidx;

mean_px = px(midx);
mean_gx = gx(midx);
% mean_dx = dx(midx);

% mean_px = px;
% mean_gx = gx;

err = 1e-04;

abs_gx = abs(mean_gx) + err;
abs_px = abs(mean_px) + err;

w = abs_px ./ abs_gx;
w = min(w, 1);
% w = imfilter(w, kg, 'symmetric');

pidx = mean_gx.*mean_px <= 0;
% pidx = bsxfun(@times, gx, mean_px) <= 0;
% r_gx = gx;
r_gx = bsxfun(@times, gx, w);
% r_gx(pidx) = 0;
r_gx(repmat(pidx, [1 1 nch])) = 0; 

if do_normalize
  gx0 = imfilter(I0, [0 -1 1], 'replicate');
  p = gaussian1d(gx(midx), ss);
  q = gaussian1d(gx0(midx), ss);

  a = (p.*q + err) ./ (p.*p + err);
  a = min(a, 1);
  a = gaussian1d(a, ss);
  % a = imfilter(a, kg, 'symmetric');
  r_gx = bsxfun(@times, r_gx, a);
end
% wx = w;
wx = repmat(w, [1 1 nch]);

end