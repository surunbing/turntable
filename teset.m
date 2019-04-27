clc, clear
close all

K = 355.1054;
taue = 0.001668255197954;
taum = 2.196498493001917;

G_model = tf(K, [taue * taum taue + taum 1 0]);

T = 1 / (10 * 2 * pi) / 10; 
Inertial = tf(1, [T, 1]);

fre = 30 * 2 * pi;
m = sin((20) / 180 * pi);
alpha = (1 - m) / (1 + m);
tau = 1 / (sqrt(alpha) * fre);
T = alpha * tau;
G_advance1 = tf([tau, 1], [T, 1]);

fre = 30 * 2 * pi;
m = sin((15) / 180 * pi);
alpha = (1 - m) / (1 + m);
tau = 1 / (sqrt(alpha) * fre);
T = alpha * tau;
G_advance2 = tf([tau, 1], [T, 1]);

fre = 30 * 2 * pi;
m = sin((10) / 180 * pi);
alpha = (1 - m) / (1 + m);
tau = 1 / (sqrt(alpha) * fre);
T = alpha * tau;
G_advance3 = tf([tau, 1], [T, 1]);

alpha = 3;
fre = 0.2 * 2 * pi;
tau = 1 / (sqrt(alpha * fre));
G_later1 = tf([tau, 1], [alpha * tau, 1]);

alpha = 3;
fre = 0.5 * 2 * pi;
tau = 1 / (sqrt(alpha * fre));
G_later2 = tf([tau, 1], [alpha * tau, 1]);

alpha = 3;
fre = 1 * 2 * pi;
tau = 1 / (sqrt(alpha * fre));
G_later3 = tf([tau, 1], [alpha * tau, 1]);

gain = 27 * 53.0884;

G = G_model * Inertial * G_advance1 * G_advance2 * G_advance3 * G_advance1 * G_advance1 * G_later1 * G_later2 * G_later3 * gain;

margin(G);
grid on

bode(G / (1 + G));
grid on

frequence = linspace(1, 10, 10) * 2 * pi;
[mag, phi] = bode(G / (1 + G), frequence);
mag = reshape(mag, [10, 1]);
phi = reshape(phi, [10, 1]);
