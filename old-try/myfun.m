function [q] = myfun(x, data, f1, e, T, alpha, num)
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
f2 = x(1);
fre = x(2);
tau = 1 / (sqrt(alpha) * fre);
f = data.fre(num);
complex_wcp = (complex(f1 * f1 - f * f, e * T * f) / f1 / f1) / (complex(f2 * f2 - f * f, T * f) / f2 / f2) * complex(1, tau * f) / complex(1, f * alpha * tau);
q(1) = 20 * log10(abs(complex_wcp)) + data.mag(num);
q(2) = angle(complex_wcp) + data.phi(num) + 135;
end

