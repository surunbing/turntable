clc, clear
close all

nRp = 15;
nCon = 25;
nCRp = 20;
ratio = 3;


bandwidth = 17 * 2 * pi;

%% Object
K = 1.56 * 180 / pi;
taue = 0.0039035;
taum = 0.984871194396488;
G = tf(K, [taue * taum, taum, 1, 0]);
%% ������Ի���
T = 1 / (bandwidth / 2 / pi) / 15; 
Inertial = tf(1, [T, 1]);
G = G * Inertial;

%% ����Ƶ�ʵ�
fre = logspace(-1, 2.3) * 2 * pi;  %% 1 - 100
[mag, phi] = bode(G, fre);

data.fre = fre;
data.mag = zeros(length(fre), 1);
data.phi = zeros(length(fre), 1);

for i = 1:length(fre)
    data.mag(i) = 20 * log10(mag(1, 1, i));
    data.phi(i) = phi(1, 1, i);
end

%% ������� ��ͼ
% figurename('bode');
% subplot 211;
% semilogx(data.fre, data.mag, 'b-');
% grid on
% subplot 212;
% semilogx(data.fre, data.phi, 'b-');
% grid on

%% ����λ
phi_advance = 30;
phi_advance_margin = 0;
[advance, ratio] = FillPhase(data, ratio, phi_advance, phi_advance_margin, bandwidth);

G_P = advance.gain * G * advance.P;
% figurename('ǰ��ͨ��');
% margin(G_P);
% grid on;
% figurename('��ǰ����');
% bode(advance.P);
% grid on;

G_C = G_P / (1 + G_P);
figurename('�ջ�');
bode(G_C);
grid on

later = Boostedlfgain(data, advance, 0.01 * 2 * pi, 5);
G_PP = G_P * later.gain * later.G_later;
% later2 = Boostedlfgain(data, advance, 0.1 * 2 * pi, 1.1);
% G_PP = G_PP * later2.gain * later2.G_later;
% figurename('�ٺ�ͨ��');
% margin(G_PP);
% grid on;

trap = trapfilter(data, advance, later, 7 * 2 * pi, 20, 1.5);
G_PPP = G_PP * trap.G;
figurename('��ʼ����');
margin(G_PPP);
grid on

% % %% ����ָ��Ƶ�ʵĿ����������
[data, data_con] = UpdateOptdata(G_P, bandwidth, ratio, nRp, nCon, nCRp);

%% �����Ż�
series.real_pole = 0;   % ����  a
series.real_zero = 0;   % ���  a
% series.complex_pole = 1;    % ������  a+bj
% series.complex_zero = 1;    % �����  a+bj
series.trap = 1;
series.lead = 1;            % ����    a,b
series.count = series.real_pole + series.real_zero + series.trap * 4 + series.lead * 2;

% Լ��
lb = zeros(series.count + 1, 1) + 0.001;
lb(1) = 1;
lb(2) = 0.01;
lb(3) = 0.05;
lb(4) = 0.01;
lb(5) = 0.1;
lb(6) = 0.005;
lb(7) = 0.005;
% lb(8) = 0.01;
% lb(9) = 0.01;
ub = 1e19 * ones(series.count + 1, 1);
ub(2) = 500;
ub(3) = 500;
ub(4) = 1000;
ub(5) = 1000;
ub(6) = 1000;
ub(7) = 1000;
% ub(8) = 100;
% ub(9) = 100;

%% ȷ����ҪƵ��
start = [later.gain , trap.e, trap.T, trap.f, trap.f, later.alpha, later.fre];

% start = [later.gain * later2.gain, trap.e, trap.T, trap.f, trap.f, later.alpha, later.fre, later2.alpha, later2.fre];

options = optimset('Algorithm','interior-point', 'Hessian', 'bfgs', 'MaxFunEvals', 60000, 'MaxIter', 2000);
% options = optimset('Algorithm','sqp','MaxIter',1600);
tic
[X, fval, exitflag] = fmincon(@(x)MY_costfunction(x, data, series, G_P, nRp)...
    , start, [], [], [], [], lb, ub, @(x)nonlcon1(x, data_con, series, G_P, bandwidth, nCon, nCRp, 170, 6), options);
toc

P = GetTf(X, series);

figurename('������');
bode(P);
grid on

figurename('bihuan');
bode(G_P * P / (1 + G_P * P));
grid on

figurename('margin');
margin(G_P * P);
grid on



%% �ڶ������ ����һ�����ڣ���ʼ������һ������
% series1.real_pole = 0;   % ����  a
% series1.real_zero = 0;   % ���  a
% series1.complex_pole = 1;    % ������  a+bj
% series1.complex_zero = 1;    % �����  a+bj
% series1.lead = 0;            % ����    a,b
% series1.count = series.real_pole + series.real_zero + series.complex_pole * 2 + series.complex_zero * 2 + series.lead * 2;
% P = GetTf(X, series1) / X(1);
% G_P = G_P * P;
% [data, data_con] = UpdateOptdata(G_P, bandwidth, ratio, nRp, nCon, nCRp);
% 
% %% �����Ż�
% series.real_pole = 0;   % ����  a
% series.real_zero = 0;   % ���  a
% series.complex_pole = 0;    % ������  a+bj
% series.complex_zero = 0;    % �����  a+bj
% series.lead = 2;            % ����    a,b
% series.count = series.real_pole + series.real_zero + series.complex_pole * 2 + series.complex_zero * 2 + series.lead * 2;
% 
% % Լ��
% lb = zeros(series.count + 1, 1) + 0.001;
% lb(1) = 1;
% lb(2) = 0.01;
% lb(3) = 0.005;
% lb(4) = 0.01;
% lb(5) = 0.005;
% ub = 500 * ones(series.count + 1, 1);
% ub(2) = 1000;
% ub(3) = 100;
% ub(4) = 1000;
% ub(5) = 100;
% 
% %% ȷ����ҪƵ��
% start = [X(1) * later2.gain, X(6), X(7), later2.alpha, later2.fre];
% % start(1) = start(1) * later2.gain;
% 
% 
% % start = [later.gain, trap.poles(1), trap.poles(2), trap.zeros(1), trap.zeros(2), later.alpha, later.fre, later2.alpha, later2.fre];
% 
% 
% options = optimset('Algorithm','interior-point', 'Hessian', 'bfgs', 'MaxFunEvals', 60000, 'MaxIter', 2000);
% % options = optimset('Algorithm','sqp','MaxIter',1600);
% tic
% [X, fval, exitflag] = fmincon(@(x)MY_costfunction(x, data, series, G_P, nRp)...
%     , start, [], [], [], [], lb, ub, @(x)nonlcon1(x, data_con, series, G_P, bandwidth, nCon, nCRp, 175, 6.1), options);
% toc
% 
% P = GetTf(X, series);
% 
% figurename('������2');
% bode(P);
% grid on
% 
% figurename('bihuan2');
% bode(G_P * P / (1 + G_P * P));
% grid on
% 
% figurename('margin2');
% margin(G_P);
% grid on
% hold on
% margin(G_P * P);

autoArrangeFigures












