% G 对象
function [c,ceq] = nonlcon1(x, data, series, P, bandwidth, nCon, nCRp, phi_reg, mag_max)

[mag, phi] = GetMagPhi(x, series, data.fre);
mag = mag + data.mag;
phi = phi +data.phi;
nCount = nCon;
c = zeros(nCon + 7 + nCRp, 1);
for i = 1 : nCon
    c(i) = - phi(i) - phi_reg;
end
c(nCon + 1) = mag(nCon + 3) - 0.4;
c(nCon + 2) = -0.4 - mag(nCon + 3);
c(nCon + 3) = phi(nCon + 3) + 130;
c(nCon + 4) = -140 - phi(nCon + 3);
c(nCon + 5) = - mag(nCon + 1);
c(nCon + 6) = mag(nCon + 2);
% c(nCon + 7) = 1 - (mag(nCon + 1) - mag(nCon + 2)) / (log10(data.fre(nCon + 2) - data.fre(nCon + 1))); 
Mag = mag;
Phi = phi;
complex_bode = 10 .^ (Mag ./ 20) .* complex(cos(Phi ./ 180 .* pi), sin(Phi ./ 180 .* pi));
% 闭环
complex_bode = complex_bode ./ (1 + complex_bode);
%转成bode形式
Mag = log10(abs(complex_bode)) * 20;
%Mag_max = max(Mag);
for i = 1 : nCRp
    c(nCount + 6 + i) = Mag(nCon + 3 + i) - 6;
end
c(nCon + 7 + nCRp) = x(2) / x(3) - 1;
ceq = [];
