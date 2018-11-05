% G 对象
function [c,ceq] = nonlcon1(x, data, series, P, bandwidth, fre)

[mag, phi] = GetMagPhi(x, series, fre);
nCount = length(fre);
c = zeros(nCount * 2, 1);
for i = 1 : nCount - 2
    c(i) = phi(i) - 180;
end
c(nCount - 1) = -mag(nCount - 1);
c(nCount) = mag(nCount);

% 闭环增益
Mag = mag + data.mag;
Phi = phi + data.phi;
complex_bode = 10 .^ (Mag ./ 20) .* complex(cos(Phi ./ 180 .* pi), sin(Phi ./ 180 .* pi));
% 闭环
complex_bode = complex_bode ./ (1 + complex_bode);
%转成bode形式
Mag = log10(abs(complex_bode)) * 20;
% Phi = angle(complex_bode) / pi * 180;
for i = nCount + 1 : nCount * 2
    c(i) = 7 - Mag(i);
end

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
c(4) = Pm - 55;
c(5) = 150 - Wmp;
c(6) = abs(phi_min) - 170;
ceq = [ ];