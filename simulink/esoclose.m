clc, clear
close all

K = 1.56 * 180 / pi;
taue = 0.0039035;
taum = 0.984871194396488;

taum2 = taum * 1;

G =  tf(K, [taue * taum2, taue + taum2, 1, 0]);
G2 = tf([taue * taum, taue + taum, 1, 0], K);


wq = 60 * 2 * pi;
% beta4 = 0.2756;
% beta1 = 0.9528;
% beta2 = 1.4539;
% beta3 = 0.7426;

beta4 = 1;
beta1 = 4;
beta2 = 6;
beta3 = 4;
Q = tf(beta4 * (wq ^ 4), [1, beta1 * wq, beta2 * wq * wq, beta3 * wq * wq * wq, beta4 * wq * wq * wq * wq]);
figurename('ESO_filter');
bode(Q);
grid on

G_eso = (G2 * G - 1) * Q;
GESO_CLOSE = 1 / (1 + G_eso);

load('controller_10_3_nolow.mat');

figurename('close');
G_close = (P * GESO_CLOSE * G) / (1 + P * GESO_CLOSE * G);
bode(G_close);
grid on

autoArrangeFigures;
