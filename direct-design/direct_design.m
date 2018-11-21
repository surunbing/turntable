function [P, G, para] = direct_design()
%% ����Ҫ�޸�г���Լ��
%% ������������
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

%% ����˫ʮָ����Ż�
start = [0.3, 405.4322];
lb = [0.00001; 0.00001];
ub = [1; inf];
% options = optimset('Algorithm','interior-point');
% options = optimset('Algorithm','sqp');
% [X, fval, exitflag] = fmincon(@(x)GetWsCost(x, T)...
%     , start, [], [], [], [], lb, ub, @(x)nlconws(x, T, pmr, kgr, Mre, wcmax), options);


%% ���ڼ���Ƶ�ʵ��Ż�
wfr = 18 * 2 * pi; 
[X, fval, exitflag] = fmincon(@(x)GetWcCost(x, T)...
    , start, [], [], [], [], lb, ub, @(x)nlconwc(x, T, wfr, kgr, Mre, pmr), options);

%% ������λԣ�ȵ��Ż�
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
figurename('�ջ�����');
G = tf(omegan * omegan, conv([1, 2 * xi * omegan, omegan * omegan], [T, 1]));
bode(G);
grid on

%% ��������
a = omegan * omegan * conv([taue, 1], [taum, 1]);
b = K * [T, 2 * T * xi * omegan + 1, omegan * (T * omegan + 2 * xi)];
P = tf(a, b);
G = tf(K, [taum * taue, taue + taum, 1, 0]);
G_P = G * P;
figurename('��������');
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