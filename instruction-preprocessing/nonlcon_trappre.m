function [c, ceq] = nonlcon_trap(x, data, num, phi_reg, mag_reg)

%%µÃµ½
e = zeros(num, 1);
T = zeros(num, 1);
f = zeros(num, 1);
for i = 1 : num
    e(i) = abs(x(i * 3 - 2));
    T(i) = abs(x(i * 3 - 1));
    f(i) = abs(x(i * 3));
end
frequence = data.fre;
complex_trap = ones(length(frequence), 1);
for i = 1 : num
    complex_trap = complex_trap .* complex(f(i) * f(i) - frequence .* frequence, e(i) * T(i) * frequence) ./ complex(f(i) * f(i) - frequence .* frequence, T(i) * frequence);
end

complex_P = data.mag .* complex(cos(data.phi / 180 * pi), sin(data.phi / 180 * pi));
complex_P = complex_P .* complex_trap;
mag = 20 * log10(abs(complex_P));
phi = angle(complex_P) / pi * 180;

c = zeros(length(frequence), 1);

for i = 1 : length(frequence)
    c(4 * i - 3) = -phi_reg - phi(i);
    c(4 * i - 2) = phi(i) - phi_reg;
    c(4 * i - 1) = mag(i) - mag_reg;
    c(4 * i) = - 0.0 - mag(i);
end

ceq = 0;
end