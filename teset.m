clc, clear
close all

data = load('ET2050416.csv');

res.fre = data(:, 1) * 2 * pi;
res.mag = data(:, 2);
res.phi = data(:, 3);

res.complex = res.mag .* complex(cos(res.phi / 180 * pi), sin(res.phi / 180 * pi));


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
fre = 0.1 * 2 * pi;
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

gain = 50.91;

G = G_model;
P = Inertial * G_advance1 * G_advance2 * G_advance1 * G_advance1 * gain;
G_auxiliary = tf(1, conv([1 / (200 * 2 * pi), 1], [1 / (200 * 2 * pi), 1]));
forward_app = tf(K, [taue * taum taue + taum 1]);
forward = 0.3 / forward_app * G_auxiliary * tf([1, 0], [0.0075, 1]);

[magp, phip] = bode(P, res.fre);
magp = reshape(magp, [length(res.fre), 1]);
phip = reshape(phip, [length(res.fre), 1]);
complexp = magp .* complex(cos(phip / 180 * pi), sin(phip / 180 * pi));

cc = complexp .* res.complex;
cc = cc ./ (1 + cc);

magc = abs(cc);
phic = angle(cc)* 180 / pi;



% figure(1);
% margin(P * G);
% grid on
% 
% figure(2);
% bode(P * G / (1 + P * G));
% grid on
% 
% figure(3);
% bode((P * G + G * forward)/ (1 + P * G));
frequence = linspace(1, 15, 15) * 2 * pi;
[mag, phi] = bode((P * G + G * forward)/ (1 + P * G), frequence);
% [mag, phi] = bode((P * G)/(1 + P * G), frequence);
mag = reshape(mag, [15, 1]);
phi = reshape(phi, [15, 1]);
phi = phi;






a = [1, 0];
b = [0.001, 1];
TSp = 0.0005;
[dNumf, dDenf] = c2dm(a, b, TSp, 'tustin');
    
ccc = 3.0875 * complex(cos(-137.591 / 180 * pi), sin(-37.591 / 180 * pi));
ccc = ccc / (1 + ccc);
abs(ccc)
angle(ccc) / pi * 180
