function [traX,traY,traZ] = gen_tubulin3d(walkLength, dk_lat, dk_axi, k0_lat, k0_axi, ...
    kmin_lat, kmax_lat, kmin_axi, kmax_axi, v0, x0)

W_lat = dk_lat.* randn(walkLength, 1);
Wn_lat = zeros(walkLength, 1);
for n = 1:walkLength
    Wn_lat(n) = sum(W_lat(1:n));
end

W_axi = dk_axi.* randn(walkLength, 1);
Wn_axi = zeros(walkLength, 1);
for n = 1:walkLength
    Wn_axi(n) = sum(W_axi(1:n));
end

kn_lat = exp(2*pi*1j* clip(k0_lat + Wn_lat, kmin_lat, kmax_lat));
kn_axi = exp(2*pi*1j* clip(k0_axi + Wn_axi, kmin_axi, kmax_axi));

vn = zeros(walkLength, 3);
for n = 1:walkLength
    phi = angle(sum(kn_lat(1:n)));
    theta = angle(sum(kn_axi(1:n)));    
    vn(n,:) = v0 .* [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
end

xn = zeros(walkLength, 3);
for n = 1:walkLength
    xn(n,:) = x0 + sum(vn(1:n,:), 1);
end

traX = [x0(1);xn(:,1)];
traY = [x0(2);xn(:,2)];
traZ = [x0(3);xn(:,3)];

end
