clc, clear
close all

%% 获取最初设计的内容
[P, G, para] = direct_design();

%% 波形修正，加入滞后环节降低剪切频率
%% 检查波形，加入陷波滤波器提高增益
%% 加入迟后环节可以引起频带以内较高的部分有一定的谐振，闭环时由陷波滤波器引起的幅频损失可以到允许范围内
%% 不引起低频增益的损失

% 首先添加之后环节, 检查此wc是否可以满足此种方法
wc_r = 250;
[mag, phi] = bode(P * G, wc_r);
wc_max = fsolve(@(x)myfun(x, para.xi, para.omegan, para.T), 100);
if wc_max > wc_r
%     wc_r = wc_r;
else
    wc_r = wc_max;
end

phi_n_min = -130 - phi;
phi_n_max = -128 - phi;
tic
% 得到迟后环节计算 
start = [wc_r / 2, 3];
lb = [1; 5];
ub = [wc_r; 15];
% options = optimset('Algorithm','interior-point');
options = optimset('Algorithm','sqp');
[X, fval, exitflag] = fmincon(@(x)GetAlphacost(x)...
    , start, [], [], [], [], lb, ub, @(x)nonlcon(x, wc_r, phi_n_min, phi_n_max), options);
toc

fre = X(1);
alpha = X(2);
tau = 1 / (sqrt(alpha) * fre);
T = alpha * tau;
G_later = tf([tau, 1], [T, 1]);

%优化增益K，使得剪切频率
count = 50;
K = linspace(1, 5, count);
Wc = zeros(count, 1);
Pm = zeros(count, 1);
Phi = zeros(count, 1);
for i = 1 : count
    [Gm, pm, Wgm, Wpm] = margin(K(i) * G_later * P * G);
    Wc(i) = Wpm;
    Pm(i) = pm;
    if pm < 50 || pm > 55
        Pm(i) = 0;
    end
    if Wpm > wc_r || Wpm < (wc_r - 50)
        Wc(i) = 0;
    end
    [~, Phi(i)] = bode((K(i) * G_later * P * G) / (1 + K(i) * G_later * P * G), para.dt);
end

% 计算增益K


figurename('开环特性');
K = P * G * G_later * 2.3;
margin(K);
grid on
figurename('闭环特性');
bode(K / (1 + K));
grid on

autoArrangeFigures;
% close all
% p = tf([tau, 1], [(alpha + 1) * tau, 2]);
% bode(p);
% grid on


