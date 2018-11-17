function [c, ceq] = nonlcon_trap(x, wc, pm_cost)

%%µÃµ½
e = x(1);
T = x(2);
f = x(3);
e1 = x(4);
T1 = x(5);
f1 = x(6);

complex_trap = complex(f * f - wc * wc, e * T * wc) / complex(f * f - wc * wc, T * wc);
complex_trap = complex_trap * complex(f1 * f1 - wc * wc, e1 * T1 * wc) / complex(f1 * f1 - wc * wc, T1 * wc);
phi = abs(angle(complex_trap) / pi * 180);

c(1) = phi - pm_cost;
ceq = 0;

end

