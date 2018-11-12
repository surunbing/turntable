clc, clear
close all

%% 还需要修改谐振峰约束

%% 给定对象特性

K = 434;
taum = 0.67;
taue = 0.0035;

T = 0.0014;
wfr = 10 * 2 * pi;
kgr = 5;
Mre = 6;
wcmax = 315;
start = [0.8, 405.4322];
lb = [0; 0];
ub = [1; inf];
options = optimset('Algorithm','interior-point');
% options = optimset('Algorithm','sqp');



[X, fval, exitflag] = fmincon(@(x)GetPmCost(x, T)...
    , start, [], [], [], [], lb, ub, @(x)nlcon(x, T, wfr, kgr, Mre, wcmax), options);
% [x, fval, exitflag] = fmincon(@(x)GetPmCost(x, T), [0.2807, 405.4322], [], [], [], [], [0, 0], [inf, inf], @(x)nlcon(x, T, wfr, kgr, Mre, wcmax), options);

x = X;
kg = 20 * log10(GetGm(x, T));
wc = GetWc(x, T);
pm = GetPm(x, T);
Mr = 20 * log10(GetRp(x, T));
dt = GetDt(x, T);

omegan = x(2);
xi = x(1);
G = tf(omegan * omegan, conv([1, 2 * xi * omegan, omegan * omegan], [T, 1]));
bode(G);
grid on

