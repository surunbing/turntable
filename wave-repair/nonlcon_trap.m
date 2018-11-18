function [c, ceq] = nonlcon_trap(x, wc, pm_cost, data)

%%�õ�
e = x(1);
T = x(2);
f = x(3);
e1 = x(4);
T1 = x(5);
f1 = x(6);
frequence = data.fre;

c = zeros(length(frequence) * 2 + 2, 1);

complex_trap_wc = complex(f * f - wc * wc, e * T * wc) / complex(f * f - wc * wc, T * wc);
complex_trap_wc = complex_trap_wc * complex(f1 * f1 - wc * wc, e1 * T1 * wc) / complex(f1 * f1 - wc * wc, T1 * wc);
phi = abs(angle(complex_trap_wc) / pi * 180);

c(1) = phi - pm_cost;

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
    c(2 * i) = -10 - phi(i);
    c(2 * i + 1) = phi(i) - 10;
end
c(length(frequence) * 2 + 2) = mag(length(frequence)) - 1;
    
    
    
ceq = 0;

end

