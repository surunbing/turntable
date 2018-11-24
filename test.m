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

% t = 0: 0.0005: 30;
% u = ones(length(t), 1) * 3;
% out = lsim(P * G / (1 + P * G), u, t);
% figurename('阶跃1');
% plot(t, u, 'b');
% hold on
% grid on
% plot(t, out, 'r');


%% 是否需要加入能否设计出的评估
[trap, later, bfailure, data_check] = wave_repair(P, G, para, 270, 0, bandwidth);
K = P * G * later.G;
for i = 1 : trap.num
    K = K * trap.G(i);
end
figurename('陷波滤波器');
margin(K);
grid on
figurename('陷波滤波器闭环');
bode(K / (1 + K));
grid on

% t = 0: 0.0005: 30;
% u = ones(length(t), 1) * 1;
% out = lsim(K / (1 + K), u, t);
% figurename('阶跃');
% plot(t, u, 'b');
% hold on
% grid on
% plot(t, out, 'r');

Design_Lowgain

figurename('低频增益');
K = K * Glow;
margin(K);
grid on
figurename('低频闭环');
bode(K / (1 + K));
grid on


autoArrangeFigures;
