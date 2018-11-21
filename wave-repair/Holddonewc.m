function [later, fval, exitflag] = Holddonewc(P, G, para, data, wc_r, phi_dist, option)

%% 获取最初设计的内容

%% 波形修正，加入滞后环节降低剪切频率
%% 检查波形，加入陷波滤波器提高增益
%% 加入迟后环节可以引起频带以内较高的部分有一定的谐振，闭环时由陷波滤波器引起的幅频损失可以到允许范围内
%% 不引起低频增益的损失

% 首先添加之后环节, 检查此wc是否可以满足此种方法   wc的选择位置讨论，如果期望相角裕度 
% wc_r = wc_up - 5;
% phi_dist = 122;
% wc_max = fsolve(@(x)myfun(x, para.xi, para.omegan, para.T, phi_dist), 100);
% if wc_max > wc_r
% %     wc_r = wc_r;
% else
%     wc_r = wc_max;
% end
[~, phi] = bode(P * G, wc_r);
phi_n_min = -phi_dist - phi;
phi_n_max = -phi_dist + 3 - phi;

%%计算无滞后的相频特性
frequence = [wc_r, wc_r - 50, data];%para.dt];
[mag, phi] = bode(P * G, frequence);
c_data = zeros(length(frequence), 1);
for i = 1 : length(frequence)
    c_data(i) = mag(1, 1, i) * complex(cos(phi(1, 1, i) / 180 * pi), sin(phi(1, 1, i) / 180 * pi));
end

tic
% 得到迟后环节计算 
start = [wc_r / 2, 6, 2];
lb = [1; 1.5; 0.5];
ub = [wc_r; 15; 5];
% options = optimset('Algorithm','interior-point');
options = optimset('Algorithm','sqp');
[X, fval, exitflag] = fmincon(@(x)GetAlphacost(x, c_data, frequence)...
    , start, [], [], [], [], lb, ub, @(x)nonlcon(x, wc_r, phi_n_min, phi_n_max, c_data, frequence), options);
toc

% [c, ceq] = nonlcon(X, wc_r - 20, phi_n_min, phi_n_max, c_data, frequence);

fre = X(1);
alpha = X(2);

tau = 1 / (sqrt(alpha) * fre);
T = alpha * tau;
later.G = tf([tau, 1], [T, 1]) * X(3);
later.K = X(3);
later.fre = X(1);
later.alpha = X(2);
% P = P * G_later * X(3);
% figurename('陷波前');
% margin(P * G);
% grid on
% trap = trapdesign(P, G, 18 * 2 * pi);
% figurename('陷波前闭环');
% bode(P * G / (1 + P * G));
% grid on




% K = 1.56 * 180 / pi;
% taue = 0.0039035;
% taum = 0.984871194396488 * 15;
% G = tf(K, [taum * taue, taue + taum, 1, 0]);


% figurename('开环特性');
% K = P * trap.G1 * G * trap.G2 * trap.G3;
% margin(K);
% grid on
% figurename('闭环特性');
% bode(K / (1 + K));
% grid on




% %% 增加陷波环节
% e = 1.5;
% T = 15;
% f1 = 70;
% 
% trap = tf([1, e * T, f1 * f1], [1, T, f1 * f1]);
% K = K * trap;
% 
% e = 2.0;
% T = 15;
% f1 = 88;
% trap = tf([1, e * T, f1 * f1], [1, T, f1 * f1]);
% K = K * trap;
% 
% e = 1.8;
% T = 15;
% f1 = 113;
% trap = tf([1, e * T, f1 * f1], [1, T, f1 * f1]);
% 
% 
% 
% %% 迟后环节，增加开环增益
% alpha = 2;
% fre = 0.5;
% tau = 1 / sqrt(alpha) / fre;
% T = alpha * tau;
% G_later = tf([tau, 1], [T, 1]) * alpha;
% 
% 
% K = K * trap * G_later * G_later;
% 
% margin(K);
% grid on
% figurename('闭环特性');
% bode(K / (1 + K));
% grid on

autoArrangeFigures;
% close all
% p = tf([tau, 1], [(alpha + 1) * tau, 2]);
% bode(p);
% grid on


