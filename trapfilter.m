function [trap] = trapfilter(data, advance, later, fre, width, depth)
e = depth;
T = width;
f = fre;

a = [1, e * T, f^2];
b = [1, T, f^2];

G = tf(a, b);
trap.e = e;
trap.T = T;
trap.f = f;
trap.G = G;
zero_points = roots(a);
poles_points = roots(b);
trap.zeros = [-real(zero_points(1)), imag(zero_points(1))];
trap.poles = [-real(poles_points(1)), imag(poles_points(1))];