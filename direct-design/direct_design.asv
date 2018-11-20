function [P, G, para] = direct_design()
%% 还需要修改谐振峰约束
%% 给定对象特性
% K = 434;
% taum = 0.67;
% taue = 0.0035;

K = 1.56 * 180 / pi;
taue = 0.0039035;
taum = 0.984871194396488;

T = 0.0014 / 2;
kgr = 5;
Mre = 6;
wcmax = 550;
pmr = 35;

%% 基于双十指标的优化
start = [0.3, 405.4322];
lb = [0.00001; 0.00001];
ub = [1; inf];
% options = optimset('Algorithm','interior-point');
% options = optimset('Algorithm','sqp');
% [X, fval, exitflag] = fmincon(@(x)GetWsCost(x, T)...
%     , start, [], [], [], [], lb, ub, @(x)nlconws(x, T, pmr, kgr, Mre, wcmax), options);


%% 基于剪切频率的优化
wfr = 18 * 2 * pi; 
[X, fval, exitflag] = fmincon(@(x)GetWcCost(x, T)...
    , start, [], [], [], [], lb, ub, @(x)nlconwc(x, T, wfr, kgr, Mre, pmr), options);

%% 基于相位裕度的优化
% wfr = abs(fval); 
% % wcmax = 350;
% start = [0.5, 405.4322];
% lb = [0.00001; 0.00001];
% ub = [inf; inf];
% options = optimset('Algorithm','interior-point');
% % options = optimset('Algorithm','sqp');
% [X, fval, exitflag] = fmincon(@(x)GetPmCost(x, T)...
%     , start, [], [], [], [], lb, ub, @(x)nlconpm(x, T, wfr, kgr, Mre, wcmax), options);
% [c, ceq] = nlconpm(X, T, wfr, kgr, Mre, wcmax);
x = X;
kg = 20 * log10(GetGm(x, T));
wc = GetWc(x, T);
pm = GetPm(x, T);
Mr = 20 * log10(GetRp(x, T));
dt = GetDt(x, T);

omegan = x(2);
xi = x(1);
figurename('闭环对象');
G = tf(omegan * omegan, conv([1, 2 * xi * omegan, omegan * omegan], [T, 1]));
bode(G);
grid on

%% 开环对象
a = omegan * omegan * conv([taue, 1], [taum, 1]);
b = K * [T, 2 * T * xi * omegan + 1, omegan * (T * omegan + 2 * xi)];
P = tf(a, b);
G = tf(K, [taum * taue, taue + taum, 1, 0]);
G_P = G * P;
figurename('开环对象');
margin(G_P);
grid on

para.kg = kg;
para.dt = dt;
para.wc = wc;
para.pm = pm;
para.mr = Mr;
para.xi = x(1);
para.omegan = x(2);
para.T = T;