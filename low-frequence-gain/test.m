clc, clear
close all

load('control.mat');

figurename('��������');
margin(K);
grid on




autoArrangefigures;

 tau = 1 / (sqrt(0.7) * 2 * pi);
 GG = tf([tau, 1], [0.7 * tau, 1]);
 bode(GG);
 grid on

