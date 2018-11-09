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
taum = 0.984871194396488 * 3;
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

G_P = G * advance.P;

figurename('ǰ��ͨ��');
margin(G_P);
grid on;
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



autoArrangeFigures












