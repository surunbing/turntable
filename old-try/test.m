clc, clear
close all

e = 0.9;
T = 700;
f1 = 113;
f2 = 50;

G = tf([1, e * T, f1 * f1], [1, T, f1 * f1]);
bode(G);
grid on
hold on
G = tf([1, e * T, f1 * f1], [1, T, f2 * f2]);
bode(G);


% fre = 0.01;
% alpha = 5;
% tau = 1 / (sqrt(alpha) * fre);
% T = alpha * tau;
% G_later = tf([tau, 1], [T, 1]);
% figurename('迟后环节');
% bode(G_later);
% grid on
% 
% K = 434 / 1.2;
% taum = 0.67;
% taue = 0.0035;
% G = tf(K, [taum * taue taum + taue 1 0]);
% K = tf(164375.2706 * conv([0.0035 1], [0.67, 1]), 434 * [0.0014, 1.3187, 457.7607]);
% G_later = G_later * 3;
% figurename('开环');
% margin(G * K * G_later);
% grid on
% 
% figurename('闭环');
% bode(G * K * G_later/ (1 + G * K * G_later));
% grid on