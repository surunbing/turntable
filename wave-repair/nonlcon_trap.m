function [c, ceq] = nonlcon_trap(x, wc, pm_cost, data, num)

%%得到
e = zeros(num, 1);
T = zeros(num, 1);
f = zeros(num, 1);
for i = 1 : num
    e(i) = x(i * 3 - 2);
    T(i) = x(i * 3 - 1);
    f(i) = x(i * 3);
end
frequence = data.fre;
frequence(length(frequence)) = frequence(length(frequence)) + pi;
complex_trap = ones(length(frequence), 1);
complex_trap_wc = 1;
for i = 1 : num
    complex_trap = complex_trap .* complex(f(i) * f(i) - frequence .* frequence, e(i) * T(i) * frequence) ./ complex(f(i) * f(i) - frequence .* frequence, T(i) * frequence);
    complex_trap_wc = complex_trap_wc * complex(f(i) * f(i) - wc * wc, e(i) * T(i) * wc) / complex(f(i) * f(i) - wc * wc, T(i) * wc);
end

c = zeros(length(frequence) * 2 + 2 + num, 1);

phi = abs(angle(complex_trap_wc) / pi * 180);

c(1) = phi - pm_cost;

complex_P = data.mag .* complex(cos(data.phi / 180 * pi), sin(data.phi / 180 * pi));
complex_P = complex_P .* complex_trap;
complex_c = complex_P ./ (1 + complex_P);
mag = 20 * log10(abs(complex_c));
phi = angle(complex_c) / pi * 180;

for i = 1 : length(frequence)
    c(2 * i) = -10 - phi(i);
    c(2 * i + 1) = phi(i) - 10;
end
c(length(frequence) * 2 + 2) = mag(length(frequence)) - 0.8;

%% 求取xian滤波器的最小值
x = ((T.^2.*e)/2 + f.^2 + (T.*(e.*(e.*T.^2 + 4.*f.^2)).^(1/2))./2).^(1/2);
complex_phi = complex(f .* f - x .* x, e .* T .* x) ./ complex(f .* f - x .* x, T .* x);
phi = angle(complex_phi) / pi * 180;
c(length(frequence) * 2 + 3 : length(frequence) * 2 + 2 + num) = -20 - phi;
ceq = 0;

end

