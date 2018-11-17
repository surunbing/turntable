function [cost] = GetTrapcost(x, num, data, wc, pm_min)
cost = 0;
%	最小化代价函数和
e = x(1);
T = x(2);
f = x(3);
e1 = x(4);
T1 = x(5);
f1 = x(6);

frequence = data.fre;
complex_trap = zeros(length(frequence), 1);
for i = 1 : 1 : length(frequence)
    complex_trap(i) = complex(f * f - frequence(i) * frequence(i), e * T * frequence(i)) / complex(f * f - frequence(i) * frequence(i), T * frequence(i));
    complex_trap(i) = complex_trap(i) * complex(f1 * f1 - frequence(i) * frequence(i), e1 * T1 * frequence(i)) / complex(f1 * f1 - frequence(i) * frequence(i), T1 * frequence(i));
end
complex_P = data.mag .* complex(cos(data.phi / 180 * pi), sin(data.phi / 180 * pi));
complex_P = complex_P .* complex_trap;
complex_c = complex_P ./ (1 + complex_P);
mag = 20 * log10(abs(complex_c));
phi = angle(complex_c) / pi * 180;

for i = 1 : length(frequence)
   cost = cost + mag(i) * mag(i) * 10 + phi(i) * phi(i); 
end

end

