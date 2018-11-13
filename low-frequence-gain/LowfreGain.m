function [q] = LowfreGain()
%   增加低频增益
close all
fre = 0.5 * 2 * pi;
alpha = 2;
tau = 1 / (sqrt(alpha) * fre);
T = tau * alpha;
G = tf([tau, 1], [T, 1]);
fre = 0.3 * 2 * pi;
alpha = 2;
tau = 1 / (sqrt(alpha) * fre);
T = tau * alpha;
G1 = tf([tau, 1], [T, 1]);
bode(G);
grid on
[mag, phi] = bode(G * G1 * G1, 100);
db = 20 * log10(2);
end

