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

%% 是否需要加入能否设计出的评估
bandwidth1 = max([bandwidth +pi, para.dt]);
[later, fval, exitflag] = Holddonewc(P, G, para, bandwidth1, bandwidth * 2.2, 0);
figurename('迟后');
margin(P * G * later.G);
grid on
figurename('迟后闭环');
bode(P * G * later.G / (1 + P * G * later.G));
grid on

autoArrangeFigures;
