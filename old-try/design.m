clc, clear
close all

nRp = 15;
nCon = 25;
nCRp = 20;
ratio = 2.2;

bandwidth = 18 * 2 * pi;

%% Object
K = 45;%1.56 * 180 / pi;
taue = 0.00396;% 0.0039035;
taum = 0.09947;%0.984871194396488;
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
advance = FillPhase(data, ratio, phi_advance, phi_advance_margin, bandwidth);

G_P = G * advance.P;

figurename('ǰ��ͨ��');
margin(G_P);
grid on;
figurename('bihuan1');
bode(G_P / (1 + G_P));
grid on
% figurename('��ǰ����');
% bode(advance.P);
% grid on;
% 
% G_C = G_P / (1 + G_P);
% figurename('�ջ�');
% bode(G_C);
% grid on
% 
% later = Boostedlfgain(data, advance, 0.01 * 2 * pi, 5);
% G_PP = G_P * later.gain * later.G_later;
% % later2 = Boostedlfgain(data, advance, 0.1 * 2 * pi, 1.1);
% % G_PP = G_PP * later2.gain * later2.G_later;
% % figurename('�ٺ�ͨ��');
% % margin(G_PP);
% % grid on;
% 
% trap = trapfilter(data, advance, later, 7 * 2 * pi, 20, 1.5);
% G_PPP = G_PP * trap.G;
% figurename('��ʼ����');
% margin(G_PPP);
% grid on

[Mag, Phi] = bode(G_P, data.fre);
mag = zeros(length(data.fre), 1);
phi = zeros(length(data.fre), 1);
for i = 1 : length(data.fre)
    mag(i) = 20 * log10(Mag(1, 1, i));
    phi(i) = Phi(1, 1, i);
end
data_cur.fre = data.fre;
data_cur.mag = mag;
data_cur.phi = phi;
num = advance.num;

e = 0.566315961;
T = 500;%499.99775130;
f1 = bandwidth * ratio;
alpha = 0.0402494756585578;
f2 = 92.53771513;
fre = 18.389328242984;
tau = 1 / (sqrt(alpha) * fre);
P = tf([1, e * T, f1 * f1] / f1 / f1, [1, T, f2 * f2] / f2 / f2) * tf([tau, 1], [alpha * tau, 1]);
h1 = figurename('bode');
margin(G_P * P);
grid on
[Gm, Pm, Wcm, Wcp] = margin(G_P * P);

count = round(bandwidth / 2 / pi);
Fre = linspace(1, count, count)' * 2 * pi;
[mag, phi, ~] = bode(G_P, Fre);

data.fre = Fre;
data.mag = zeros(length(Fre), 1);
data.phi = zeros(length(Fre), 1);

for i = 1:length(Fre)
    data.mag(i) = 20 * log10(mag(1, 1, i));
    data.phi(i) = phi(1, 1, i);
end
x = [e, T, f1, f2, alpha, fre];
cost = MY_costfunction(x, data);

h2 = figurename('bihuan');
bode(G_P * P / (1 + G_P * P));
grid on
N = 20;
e = linspace(0.01, 1, N);
alpha = logspace(log10(0.005), log10(0.15), N);
res = zeros(N, N, 2);
Cost = zeros(N, N);
for i = 1 : N
    for j = 1 : N
%         p = myfun([f2, fre], data_cur, f1, e, T, alpha, num);
        res(i, j, :) = fsolve(@(x)myfun(x, data_cur, f1, e(i), T, alpha(j), num), [f1 * 1.5, bandwidth / 2 / pi]');
        f2 = res(i, j, 1);
        fre = res(i, j, 2);
%         tau = 1 / (sqrt(alpha(j)) * fre);
        x = [e(i), T, f1, f2, alpha(j), fre];
        Cost(i, j) = MY_costfunction(x, data);
        f2 = res(i, j, 1);
        fre = res(i, j, 2);
        tau = 1 / (sqrt(alpha(j)) * fre);
        P = tf([1, e(i) * T, f1 * f1] / f1 / f1, [1, T, f2 * f2] / f2 / f2) * tf([tau, 1], [alpha(j) * tau, 1]);
        [Gm, Pm, Wcm, Wcp] = margin(G_P * P);
        if Pm < 30 || Gm < 1
             Cost(i, j) = 0;
        end
%         P = tf([1, e(i) * T, f1 * f1] / f1 / f1, [1, T, f2 * f2] / f2 / f2) * tf([tau, 1], [alpha(j) * tau, 1]);
%         figure(h1);
% 
%         margin(G_P * P);
%         grid on
% 
%         figure(h2);
%         bode(G_P * P / (1 + G_P * P));
%         grid on
    end
end
i = 12;
j = 1;
f2 = res(i, j, 1);
fre = res(i, j, 2);
tau = 1 / (sqrt(alpha(j)) * fre);
x = [e(i), T, f1, f2, alpha(j), fre];
P = tf([1, e(i) * T, f1 * f1] / f1 / f1, [1, T, f2 * f2] / f2 / f2) * tf([tau, 1], [alpha(j) * tau, 1]);
Wm = 0;%tf([20, 1], [1, 20, 80089]);
figurename('margin');
margin(G_P * P);
grid on
hold on
margin(G_P * P * (1 + Wm) * 1.5);
figurename('�ջ�');
bode(G_P * P * (1 + Wm) / (1 + G_P * P * (1 + Wm)));
grid on
%% ǰ��
K = 45;%1.56 * 180 / pi;
taue = 0.00396;% 0.0039035;
taum = 0.09947;%0.984871194396488;
G_F = 1 / tf(K, [taue * taum, taum, 1, 0]) * 0.5;

K = 45;%1.56 * 180 / pi;
taue = 0.00396;% 0.0039035;
taum = 0.09947;%0.984871194396488;
G = tf(K, [taue * taum, taum, 1, 0]);
%% ������Ի���
T = 1 / (bandwidth / 2 / pi) / 15; 
Inertial = tf(1, [T, 1]);

G_C = (G_P * P + G * G_F) / (1 + G_P * P);
figurename('ǰ��');
bode(G_C);
grid on
% 
% fre = bandwidth * ratio * 2;
% alpha = 2;
% tau = 1 / (sqrt(alpha) * fre);
% T = alpha * tau;
% G_later = tf([tau, 1], [T, 1]);
% figurename('�ٺ󻷽�');
% bode(G_later);
% grid on

autoArrangeFigures


