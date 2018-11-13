clc, clear
close all

%% ����Ҫ�޸�г���Լ��

%% ������������

K = 434;
taum = 0.67;
taue = 0.0035;

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
options = optimset('Algorithm','sqp');
[X, fval, exitflag] = fmincon(@(x)GetWsCost(x, T)...
    , start, [], [], [], [], lb, ub, @(x)nlconws(x, T, pmr, kgr, Mre, wcmax), options);


%% ������λԣ�ȵ��Ż�
wfr = abs(fval); 
% wcmax = 350;
start = [0.5, 405.4322];
lb = [0.00001; 0.00001];
ub = [1; inf];
% options = optimset('Algorithm','interior-point');
options = optimset('Algorithm','sqp');
[X, fval, exitflag] = fmincon(@(x)GetPmCost(x, T)...
    , start, [], [], [], [], lb, ub, @(x)nlconpm(x, T, wfr, kgr, Mre, wcmax), options);

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

%% ��������
a = omegan * omegan * conv([taue, 1], [taum, 1]);
b = K * [T, 2 * T * xi * omegan + 1, omegan * (T * omegan + 2 * xi)];
P = tf(a, b);
G = tf(K, [taum * taue, taue + taum, 1, 0]);
G_P = G * P;
figurename('��������');
margin(G_P);
grid on


e = 2.8;
T = 11;
f1 = 113;
f2 = 113;

trap = tf([1, e * T, f1 * f1], [1, T, f1 * f1]);
% bode(G);
% grid on
% hold on
% G = tf([1, e * T, f1 * f1], [1, T, f2 * f2]);
% bode(G);


fre = 5 * 2 * pi;
alpha = 5;
tau = 1 / (sqrt(alpha) * fre);
T = alpha * tau;
G_later = tf([tau, 1], [T, 1]);
figurename('�ٺ󻷽�');
bode(G_later);
grid on

K = 434;
taum = 0.67 * 2;
taue = 0.0035;
G = tf(K, [taum * taue taum + taue 1 0]);
% K = tf(164375.2706 * conv([0.0035 1], [0.67, 1]), 434 * [0.0014, 1.3187, 457.7607]);
G_later = G_later * 1.5;
figurename('����');
margin(G * P * G_later * trap);
grid on

figurename('�ջ�');
bode(G * P * G_later * trap / (1 + G * P * G_later * trap));
grid on