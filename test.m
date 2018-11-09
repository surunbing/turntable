clc, clear
close all

e = 0.566315961;
T = 499.99775130;
f1 = 320.6384535;
f2 = 92.53771513;

G = tf([1, e * T, f1 * f1], [1, T, f1 * f1]);
bode(G);
grid on
hold on
G = tf([1, e * T, f1 * f1], [1, T, f2 * f2]);
bode(G);
