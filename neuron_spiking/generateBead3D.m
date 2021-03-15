function I = generateBead3D(dims, sigma, intensity)

%% return a 3-D image with a Gaussian-like bead at the center.
%  Params:
%   -dims: 1*3 vector, [height, width, depth] of the returned 3-D array
%   -sigma: scalar or a 2-elemnet vector indicating the xy and z sigma of the
%   Gaussian distribution.

if isscalar(sigma)
    sigma = [sigma, sigma];
end

bead = zeros(dims);
[height, width, depth] = size(bead);

center_slice = fspecial('gaussian', [height, width], sigma(1));
axial_inten = fspecial('gaussian', [1, depth], sigma(2));
axial_inten = axial_inten ./ max(axial_inten(:));

for d = 1 : depth
    tmp = center_slice;
    tmp = tmp * axial_inten(1, d);
    bead(:,:,d) = tmp;
end

bead = bead ./ max(bead(:)) * intensity;
I = bead;
