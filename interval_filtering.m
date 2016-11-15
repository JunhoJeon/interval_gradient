%   The Code is created based on the method described in the following paper:
%   [1] "Structure-Texture Decomposition of Images with Interval Gradient",
%   Hyunjoon Lee, Junho Jeon, Junho Kim, Seungyong Lee, Computer Graphics Forum (CGF), 2016.
%
%   The code and the algorithm are for non-comercial use only.
%  
%   Author: Hyunjoon Lee, Junho Jeon
%   Date  : 11/14/2015
%   Version : 1.0 
%   Copyright 2015, POSTECH

function [S] = interval_filtering(img, ss)
img = im2single(img);
% I = img; % CPU version
I = gpuArray(im2single(img)); % GPU version
% I = rgb2gray(I);

sigma_s = ss;%3;
epsi = 0.03.^2;
grad_type = 'gaussian';
Nd = 3;
N = 8;

S = I;
wx_prev = zeros(size(S), 'like', S);
wy_prev = zeros(size(S), 'like', S);
wy_prev = permute(wy_prev, [2 1 3]);

total_iter = 0;
  tic;
for ii = 0:N-1
  [gx, wx] = rescale_gradient(S, sigma_s, grad_type, false, I);
  [gy, wy] = rescale_gradient(permute(S, [2 1 3]), sigma_s, grad_type, false, permute(I, [2 1 3]));
  ngx = gx;
  ngy = gy;
  
  dwx = (wx - wx_prev).^2;
  dwy = (wy - wy_prev).^2;
  mean_dwx = mean(dwx(:));
  mean_dwy = mean(dwy(:));
%   disp([mean_dwx, mean_dwy]);
  
  if min(mean_dwx, mean_dwy) <= 0.002
    break; 
  end
  
  S = I;
  for i = 0:Nd-1
    ss_i = sigma_s * 3 * sqrt(3) * 2^(Nd - (i + 1)) / sqrt(4^Nd - 1);
    ss_i = round(ss_i);
    
    if ss_i == 0, break; end
    
    S = guided_reconstruct_1d(S, ngx, wx, ss_i, epsi);
    S = permute(S, [2 1 3]);
    S = max(0, min(S, 1));
    
    S = guided_reconstruct_1d(S, ngy, wy, ss_i, epsi);
    S = permute(S, [2 1 3]);
    S = max(0, min(S, 1));
  end
  
  wx_prev = wx;
  wy_prev = wy;
  
%   epsi = epsi / 2;
%   sigma_s = max(1, sigma_s / sqrt(2));
  
%   imwrite(gather(S), num2str(ii, 'res_iter%d.jpg'), 'quality', 85);
  figure(101), imshow(S);
  drawnow;
  
  total_iter = total_iter + 1;
  
end
toc;
figure(101),imshow(S), drawnow;
disp(['total_iter = ' num2str(total_iter)]);
end