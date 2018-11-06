% clc, clear
% close all
cost = 0;
% % 计算频点处的特性
% % 开环
frequence = data.fre;
[mag, phi] = GetMagPhi(X, series, frequence);
mag = mag + data.mag;
phi = phi + data.phi;
complex_bode = 10 .^ (mag ./ 20) .* complex(cos(phi ./ 180 .* pi), sin(phi ./ 180 .* pi));
% 闭环
complex_bode = complex_bode ./ (1 + complex_bode);
%转成bode形式
mag = log10(abs(complex_bode)) * 20;
phi = angle(complex_bode) / pi * 180;

count = length(data.fre) - nRp;
for i = 1 : count
    cost = cost + mag(i) * mag(i) * 1.5 + phi(i) * phi(i);
end
for i =  1 : nRp
    cost = cost + mag(count + i) * mag(count + i) * 1 / nRp;
end