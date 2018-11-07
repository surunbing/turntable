clc, clear
close all
depth = 0.6;
width = 50;
fre = 100;

e = depth;
T = width;
f = fre;

a = [1, e * T, f^2];
b = [1, T, f^2 * 0.5];

G = tf(a, b);
bode(G);
grid on