clc, clear;
close all

% K = 434;
% taum = 0.67;
% taue = 0.0035;

K = 1.56 * 180 / pi;
taue = 0.0039035;
taum = 0.984871194396488;

bandwidth = 18 * 2 * pi;
wc_max = 650;
[P, G, para] = direct_design(bandwidth, wc_max, K, taum, taue);

figurename('直接设计开环');
margin(P * G);
grid on

figurename('直接设计闭环');
bode(P * G / (1 + P * G));
grid on

Design_Lowgain

%% 是否需要加入能否设计出的评估
[trap, later, bfailure, data_check] = wave_repair(P * Glow, G, para, 250, 0, bandwidth);
K = P * G * later.G * Glow;
for i = 1 : trap.num
    K = K * trap.G(i);
end
figurename('陷波滤波器');
margin(K);
grid on
figurename('陷波滤波器闭环');
bode(K / (1 + K));
grid on


% Design_Lowgain
% 
% figurename('低频增益');
% K = K * Glow;
% margin(K);
% grid on
% figurename('低频闭环');
% bode(K / (1 + K));
% grid on


autoArrangeFigures;
