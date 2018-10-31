%% 输入 零点个数 极点个数 复数零点  复数极点 
function cost = MY_costfunction(x, data, series, P)
cost = 0;
count = length(data.fre);

G = GetTf(x, series);
 
[Mag, Phi] = bode(G, data.fre);
mag = zeros(count, 1);
phi = zeros(count, 1);
for i = 1 : count
    mag(i) = 20 * log10(Mag(1, 1, i));
    phi(i) = Phi(1, 1, i);
end
% 计算频点处的特性
% 开环
mag = mag + data.mag;
phi = phi + data.phi;
complex_bode = 10 .^ (mag ./ 20) .* complex(cos(phi ./ 180 .* pi), sin(phi ./ 180 .* pi));
% 闭环
complex_bode = complex_bode ./ (1 + complex_bode);
%转成bode形式
mag = log10(abs(complex_bode)) * 20;
phi = angle(complex_bode) / pi * 180;

[mag_close, phi_close] = bode(G * P / (1 + G * P));
Mag_max = max(mag_close);

for i = 1 : count
    cost = cost + mag(i) * mag(i) * 1.5 + phi(i) * phi(i);
end

cost = cost + Mag_max * Mag_max;
