clc, clear
close all

nRp = 15;
nCon = 25;
nCRp = 20;
ratio = 3;


bandwidth = 15 * 2 * pi;

%% Object
K = 1.56 * 180 / pi;
taue = 0.0039035;
taum = 0.984871194396488;
G = tf(K, [taue * taum, taum, 1, 0]);
%% 加入惯性环节
T = 1 / (bandwidth / 2 / pi) / 15; 
Inertial = tf(1, [T, 1]);
G = G * Inertial;

%% 给定频率点
fre = logspace(-1, 2.3) * 2 * pi;  %% 1 - 100
[mag, phi] = bode(G, fre);

data.fre = fre;
data.mag = zeros(length(fre), 1);
data.phi = zeros(length(fre), 1);

for i = 1:length(fre)
    data.mag(i) = 20 * log10(mag(1, 1, i));
    data.phi(i) = phi(1, 1, i);
end

%% 名义对象 绘图
% figurename('bode');
% subplot 211;
% semilogx(data.fre, data.mag, 'b-');
% grid on
% subplot 212;
% semilogx(data.fre, data.phi, 'b-');
% grid on

%% 补相位
phi_advance = 30;
phi_advance_margin = 0;
[advance, ratio] = FillPhase(data, ratio, phi_advance, phi_advance_margin, bandwidth);

G_P = advance.gain * G * advance.P;
% figurename('前向通道');
% margin(G_P);
% grid on;
% figurename('超前环节');
% bode(advance.P);
% grid on;

G_C = G_P / (1 + G_P);
figurename('闭环');
bode(G_C);
grid on

later = Boostedlfgain(data, advance, 0.01 * 2 * pi, 5);
G_PP = G_P * later.gain * later.G_later;
later2 = Boostedlfgain(data, advance, 0.01 * 2 * pi, 1.05);
% G_PP = G_PP * later2.gain * later2.G_later;
% figurename('迟后通道');
% margin(G_PP);
% grid on;

trap = trapfilter(data, advance, later, 7 * 2 * pi, 20, 1.5);
G_PPP = G_PP * trap.G;
figurename('初始对象');
margin(G_PPP);
grid on

% % %% 给定指定频率的开环相关数据
[data, data_con] = UpdateOptdata(G_P, bandwidth, ratio, nRp, nCon, nCRp);

%% 整形优化
series.real_pole = 0;   % 极点  a
series.real_zero = 0;   % 零点  a
series.complex_pole = 1;    % 复极点  a+bj
series.complex_zero = 1;    % 复领点  a+bj
series.lead = 1;            % 环节    a,b
series.count = series.real_pole + series.real_zero + series.complex_pole * 2 + series.complex_zero * 2 + series.lead * 2;

% 约束
lb = zeros(series.count + 1, 1) + 0.001;
lb(1) = 1;
lb(6) = 0.01;
lb(7) = 0.005;
% lb(8) = 0.01;
% lb(9) = 0.005;
ub = 1e19 * ones(series.count + 1, 1);
ub(6) = 1000;
ub(7) = 100;
% ub(8) = 1000;
% ub(9) = 100;

%% 确定需要频点
start = [later.gain, trap.poles(1), trap.poles(2), trap.zeros(1), trap.zeros(2), later.alpha, later.fre];


% start = [later.gain, trap.poles(1), trap.poles(2), trap.zeros(1), trap.zeros(2), later.alpha, later.fre, later2.alpha, later2.fre];


options = optimset('Algorithm','interior-point', 'Hessian', 'bfgs', 'MaxFunEvals', 60000, 'MaxIter', 2000);
% options = optimset('Algorithm','sqp','MaxIter',1600);
tic
[X, fval, exitflag] = fmincon(@(x)MY_costfunction(x, data, series, G_P, nRp)...
    , start, [], [], [], [], lb, ub, @(x)nonlcon1(x, data_con, series, G_P, bandwidth, nCon, nCRp, 170, 6), options);
toc

P = GetTf(X, series);

figurename('控制器');
bode(P);
grid on

figurename('bihuan');
bode(G_P * P / (1 + G_P * P));
grid on

figurename('margin');
margin(G_P);
grid on
hold on
margin(G_P * P * later2.gain * later2.G_later);
hold on
margin(G_P * P);


%% 第二次设计 加入一个环节，起始点由上一个决定
%% 整形优化
series.real_pole = 0;   % 极点  a
series.real_zero = 0;   % 零点  a
series.complex_pole = 1;    % 复极点  a+bj
series.complex_zero = 1;    % 复领点  a+bj
series.lead = 2;            % 环节    a,b
series.count = series.real_pole + series.real_zero + series.complex_pole * 2 + series.complex_zero * 2 + series.lead * 2;

% 约束
lb = zeros(series.count + 1, 1) + 0.001;
lb(1) = 1;
% lb(2) = 50;
% lb(3) = 0.001;
% lb(4) = 10;
% lb(5) = 50;
lb(6) = 0.01;
lb(7) = 0.005;
lb(8) = 0.01;
lb(9) = 0.005;
ub = 500 * ones(series.count + 1, 1);
% ub(2) = 500;
% ub(3) = 1;
% ub(4) = 200;
% ub(5) = 500;
ub(6) = 1000;
ub(7) = 100;
ub(8) = 1000;
ub(9) = 100;

%% 确定需要频点
start = [X, later2.alpha, later2.fre];
start(1) = start(1) * later2.gain;


% start = [later.gain, trap.poles(1), trap.poles(2), trap.zeros(1), trap.zeros(2), later.alpha, later.fre, later2.alpha, later2.fre];


options = optimset('Algorithm','interior-point', 'Hessian', 'bfgs', 'MaxFunEvals', 60000, 'MaxIter', 2000);
% options = optimset('Algorithm','sqp','MaxIter',1600);
tic
[X, fval, exitflag] = fmincon(@(x)MY_costfunction(x, data, series, G_P, nRp)...
    , start, [], [], [], [], lb, ub, @(x)nonlcon1(x, data_con, series, G_P, bandwidth, nCon, nCRp, 175, 6.1), options);
toc

P = GetTf(X, series);

figurename('控制器2');
bode(P);
grid on

figurename('bihuan2');
bode(G_P * P / (1 + G_P * P));
grid on

figurename('margin2');
margin(G_P);
grid on
hold on
margin(G_P * P);

autoArrangeFigures












