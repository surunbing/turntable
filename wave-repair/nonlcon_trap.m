function [c, ceq] = nonlcon_trap(x, wc, pm_cost)

%%µÃµ½
e = x(1);
T = x(2);
f = x(3);

complex_trap = complex(f * f - wc * wc, e * T * wc) / complex(f * f - wc * wc, T * wc);
phi = abs(angle(complex_trap) / pi * 180);

c(1) = phi - pm_cost;
ceq = 0;

end

