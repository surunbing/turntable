function [c, ceq] = nonlcon(x, wc, phi_n_min, phi_n_max)

%%µÃµ½
fre = x(1);
alpha = x(2);
tau = 1 / (sqrt(alpha) * fre);
w_up = 1 / tau;
c(1) = w_up - wc;
complex_wc = complex(1, tau * wc) / complex(1, alpha * tau * wc);
phi = angle(complex_wc) / pi * 180;
c(2) = phi_n_min - phi;
c(3) = phi - phi_n_max;
ceq = 0;

end

