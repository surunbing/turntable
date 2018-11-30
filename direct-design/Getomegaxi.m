function [P, G, para, exitflag] = Getomegaxi(bandwidth, num1, T, option)
%% 还需要修改谐振峰约束
%% 给定对象特性
global parameter

kgr = parameter.kgr;
Mre = parameter.Mre;
pmr = parameter.pmr;

K = parameter.K;
taum = parameter.taum;
taue = parameter.taue;

if strcmp(option.type, 'double-ten')
    %% 基于双十指标的优化
    wcmax = num1;
    start = [0.3, 405.4322];
    lb = [0.00001; 0.00001];
    ub = [parameter.dampmax; inf];
%     options = optimset('Algorithm','interior-point');
    options = optimset('Algorithm','sqp');
    [X, fval, exitflag] = fmincon(@(x)GetWsCost(x, T)...
        , start, [], [], [], [], lb, ub, @(x)nlconws(x, T, pmr, kgr, Mre, wcmax), options);
elseif strcmp(option.type, 'wc')
    %% 基于剪切频率的优化
    wfr = bandwidth + 3 * pi; 
    start = [0.3, 405.4322];
    lb = [0.00001; 0.00001];
    ub = [parameter.dampmax; inf];
    options = optimset('Algorithm','sqp');
    [X, fval, exitflag] = fmincon(@(x)GetWcCost(x, T)...
        , start, [], [], [], [], lb, ub, @(x)nlconwc(x, T, wfr, kgr, Mre, pmr), options);
elseif strcmp(option.type, 'margin-phase')
    %% 基于相位裕度的优化
    wfr = bandwidth + pi;
    wcmax = num1;
    start = [0.3, 405.4322];
    lb = [0.00001; 0.00001];
    ub = [parameter.dampmax; inf];
%     options = optimset('Algorithm','interior-point');
    options = optimset('Algorithm','sqp');
    [X, fval, exitflag] = fmincon(@(x)GetPmCost(x, T)...
        , start, [], [], [], [], lb, ub, @(x)nlconpm(x, T, wfr, kgr, Mre, wcmax), options);
    % [c, ceq] = nlconpm(X, T, wfr, kgr, Mre, wcmax);
end
x = X;
kg = 20 * log10(GetGm(x, T));
wc = GetWc(x, T);
pm = GetPm(x, T);
Mr = 20 * log10(GetRp(x, T));
dt = GetDt(x, T);

para.kg = kg;
para.dt = dt;
para.wc = wc;
para.pm = pm;
para.mr = Mr;
para.xi = x(1);
para.omegan = x(2);
para.T = T;

omegan = x(2);
xi = x(1);

%% 开环对象
a = omegan * omegan * conv([taue, 1], [taum, 1]);
b = K * [T, 2 * T * xi * omegan + 1, omegan * (T * omegan + 2 * xi)];
P = tf(a, b);
G = tf(K, [taum * taue, taue + taum, 1, 0]);

