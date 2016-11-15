function rI = boxfilter1d(I, fr)

pI = padarray(double(I), [0 fr], 'symmetric');
cI = cumsum(pI, 2);
cI = padarray(cI, [0 1], 0, 'pre');

[~, w, ~] = size(I);

l = 1;
r = 2*fr+2;

k = (2*fr+1);

rI = cI(:, r:r+w-1, :) - cI(:, l:l+w-1, :);
rI = rI / k;
rI = cast(rI, 'like', I);

end