%% 输入 零点个数 极点个数 复数零点  复数极点 
function cost = MY_costfunction(x, data, series, P, fre_Rp)  %% Resonance peak
cost = 0;
% % 计算频点处的特性
% % 开环
frequence = [data.fre'; fre_Rp];
[mag, phi] = GetMagPhi(x, series, frequence);
mag = mag + data.mag;
phi = phi + data.phi;
complex_bode = 10 .^ (mag ./ 20) .* complex(cos(phi ./ 180 .* pi), sin(phi ./ 180 .* pi));
% 闭环
complex_bode = complex_bode ./ (1 + complex_bode);
%转成bode形式
mag = log10(abs(complex_bode)) * 20;
phi = angle(complex_bode) / pi * 180;

count = length(data.fre);
for i = 1 : count
    cost = cost + mag(i) * mag(i) * 1.5 + phi(i) * phi(i);
end
ncount = length(fre_Rp);
for i =  1 : ncount
    cost = cost + mag(count + ncount) * mag(count + ncount) * 1 / ncount;
end