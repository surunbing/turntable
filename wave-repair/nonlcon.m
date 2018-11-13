function [c, ceq] = nonlcon(x, wc, phi_n_min, phi_n_max, c_data, frequence)

%%µÃµ½
fre = x(1);
alpha = x(2);
K = x(3);
tau = 1 / (sqrt(alpha) * fre);
w_up = 1 / tau;
c(1) = w_up - wc;
complex_wc = complex(1, tau * wc) / complex(1, alpha * tau * wc);
phi = angle(complex_wc) / pi * 180;
cd_data = K * complex(1, tau * frequence) ./ complex(1, alpha * tau * frequence);
c(2) = phi_n_min - phi;
c(3) = phi - phi_n_max;
c(4) = 20 * log10(abs(cd_data(1) * c_data(1)));
c(5) = - 20 * log10(abs(cd_data(2) * c_data(2)));
ceq = 0;

end

