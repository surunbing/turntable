function [cost] = GetTrapcost(x, num, data, wc, pm_min)
cost = 0;
%	最小化代价函数和非线性环节代价函数
e = zeros(num, 1);
T = zeros(num, 1);
f = zeros(num, 1);
for i = 1 : num
    e(i) = x(i * 3 - 2);
    T(i) = x(i * 3 - 1);
    f(i) = x(i * 3);
end
frequence = data.fre;
complex_trap = ones(length(frequence), 1);
for i = 1 : num
    complex_trap = complex_trap .* complex(f(i) * f(i) - frequence .* frequence, e(i) * T(i) * frequence) ./ complex(f(i) * f(i) - frequence .* frequence, T(i) * frequence);
end
complex_P = data.mag .* complex(cos(data.phi / 180 * pi), sin(data.phi / 180 * pi));
complex_P = complex_P .* complex_trap;
complex_c = complex_P ./ (1 + complex_P);
mag = 20 * log10(abs(complex_c));
phi = angle(complex_c) / pi * 180;

for i = 1 : length(frequence)
   cost = cost + mag(i) * mag(i) * 50 + phi(i) * phi(i); 
end

end

