% G 对象
function [c,ceq] = nonlcon1(x, data, series, P, bandwidth, nCon, nCRp)

[mag, phi] = GetMagPhi(x, series, data.fre);
mag = mag + data.mag;
phi = phi +data.phi;
nCount = nCon;
c = zeros(nCon + 3, 1);
for i = 1 : nCon
    c(i) = - phi(i) - 170;
end
c(nCount + 1) = -mag(nCount + 1);
c(nCount + 2) = mag(nCount + 2);

% 闭环增益
Mag = mag;
Phi = phi;
complex_bode = 10 .^ (Mag ./ 20) .* complex(cos(Phi ./ 180 .* pi), sin(Phi ./ 180 .* pi));
% 闭环
complex_bode = complex_bode ./ (1 + complex_bode);
%转成bode形式
Mag = log10(abs(complex_bode)) * 20;
Mag_max = max(Mag);
% Phi = angle(complex_bode) / pi * 180;
c(nCount + 3) = Mag_max - 7;
ceq = [];
% 
% G = GetTf(x, series);
% [mag, phi, w] = bode(G * P);
% [Gm, Pm, Wcg, Wmp] = margin(G * P);
% wmin = abs(w - Wmp);
% num = find(wmin == min(wmin));
% phi = phi(1:num);
% phi_min = min(phi);
% c(1) = Wmp - 320;
% c(2) = 35 - Pm;
% c(3) = 1 - Gm;
% c(4) = Pm - 55;
% c(5) = 150 - Wmp;
% c(6) = abs(phi_min) - 170;
% ceq = [ ];