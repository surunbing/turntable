clc, clear
close all

bandwidth = 10 * 2 * pi;

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
figurename('bode');
subplot 211;
semilogx(data.fre, data.mag, 'b-');
grid on
subplot 212;
semilogx(data.fre, data.phi, 'b-');
grid on

%% ����λ
phi_advance = 30;
phi_advance_margin = 0;
ratio = 3;
advance = FillPhase(data, ratio, phi_advance, phi_advance_margin, bandwidth);

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

later = Boostedlfgain(data, advance, 0.01 * 2 * pi);
G_PP = G_P * later.gain * later.G_later;
% figurename('�ٺ�ͨ��');
% margin(G_PP);
% grid on;

trap = trapfilter(data, advance, later, 7 * 2 * pi, 20, 1.5);
G_PPP = G_PP * trap.G;
figurename('�ݲ�ͨ��');
margin(G_PPP);
grid on;
% figurename('�ݲ��ջ�');
% bode(G_PPP / (1 + G_PPP));
% grid on
% 
% 
% % %% ����ָ��Ƶ�ʵĿ����������
count = round(bandwidth / 2 / pi);
fre = linspace(1, count, count) * 2 * pi;
[mag, phi, w] = bode(G_P, fre);

data.fre = fre;
data.mag = zeros(length(fre), 1);
data.phi = zeros(length(fre), 1);
data.phi_rad = zeros(length(fre), 1);  
data.complex = zeros(length(fre), 1);

for i = 1:length(fre)
    data.mag(i) = 20 * log10(mag(1, 1, i));
    data.phi(i) = phi(1, 1, i);
    data.phi_rad(i) = phi(1, 1, i) / 180 * pi;
    data.complex(i) = 10 ^ (mag(1, 1, i) / 20) * complex(cos(data.phi_rad(i)), sin(data.phi_rad(i)));
end

%% �����Ż�

series.real_pole = 0;   % ����  a
series.real_zero = 0;   % ���  a
series.complex_pole = 1;    % ������  a+bj
series.complex_zero = 1;    % �����  a+bj
series.lead = 1;            % ����    a,b
series.count = series.real_pole + series.real_zero + series.complex_pole * 2 + series.complex_zero * 2 + series.lead * 2;

% Լ��
lb = zeros(series.count + 1, 1) + 0.000001;
lb(6) = 0.01;
lb(7) = 0.001;
ub = 1e19 * ones(series.count + 1, 1);
ub(6) = 1000;
ub(7) = 1000;

% start = zeros(series.count + 1, 1);
% start(1) = 1;
start = [later.gain, trap.poles(1), trap.poles(2), trap.zeros(1), trap.zeros(2), later.alpha, later.fre];
%start = [later.gain, 2, 42, -2, 41, later.alpha, later.fre];

% P = GetTf(start, series);

options = optimset('Algorithm','interior-point');
% options = optimset('Algorithm','sqp','MaxIter',1600);
tic
[X, fval, exitflag] = fmincon(@(x)MY_costfunction(x, data, series, G_P)...
    , start, [], [], [], [], lb, ub, @(x)nonlcon(x, data, series, G_P, bandwidth), options);
toc

P = GetTf(X, series);

figurename('������');
bode(P);
grid on

figurename('����');
bode(G_P * P);
grid on

figurename('bihuan');
bode(G_P * P / (1 + G_P * P));
grid on

figurename('margin');
margin(G_P);
grid on
hold on
margin(G_P * P);
autoArrangeFigures












