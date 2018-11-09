function cost = MY_costfunction(x, data)  %% Resonance peak
cost = 0;
% % 计算频点处的特性
% % 开环
nLength = length(data.fre);
mag = ones(nLength, 1);
phi = zeros(nLength, 1);
fre = data.fre;
e = x(1);
T = x(2);
f = x(3);
f1 = x(4);
res = (complex(f * f * ones(nLength, 1) - fre .* fre, e * T * fre) / f / f) ./ (complex(f1 * f1 * ones(nLength, 1) - fre .* fre, T * fre) / f1 / f1);
mag = mag .* abs(res);
phi = phi + angle(res);

alpha = x(5);
frequence = x(6);
tau = 1 / (sqrt(alpha) * frequence);
T = alpha * tau;
res = complex(ones(nLength, 1), tau * fre) ./ complex(ones(nLength, 1),  T * fre);
mag = mag .* abs(res);
phi = phi + angle(res);

mag = 20 * log10(mag);
phi = phi / pi * 180;

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