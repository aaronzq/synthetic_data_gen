function [traX,traY] = gen_tubulin(walkLength, dk, k0, kmin, kmax, v0, x0)

W = dk.* randn(walkLength, 1);
Wn = zeros(walkLength, 1);
for n = 1:walkLength
    Wn(n) = sum(W(1:n));
end

kn = exp(2*pi*1j* clip(k0 + Wn, kmin, kmax));

vn = zeros(walkLength, 1);
for n = 1:walkLength
    vn(n) = v0 * exp( 1j*angle(sum(kn(1:n))));
end

xn = zeros(walkLength, 1);
for n = 1:walkLength
    xn(n) = x0 + sum(vn(1:n));
end

traX = real([x0;xn]);
traY = imag([x0;xn]);

end

