clc, clear
syms x T e f
f(x) = (f * f - x * x) * (e - 1) * T * x / ((f * f - x * x) ^ 2 + e * T * T * x * x);
difff = diff(f(x));
res = solve(difff==0,'x')
r = res(6);
res = subs(r,[e T f],[4, 15, 100])

e = 4;
T = 15;
f = 100;
G = tf([1, e * T, f * f], [1, T, f * f]);
bode(G);
grid on;
(100^2 + (15*(16*100^2 + 3600)^(1/2))/2 + 450)^(1/2)

((T^2*e)/2 + f^2 + (T*(e*(e*T^2 + 4*f^2))^(1/2))/2)^(1/2);
((T^2*e)/2 + f^2 - (T*(e*(e*T^2 + 4*f^2))^(1/2))/2)^(1/2)