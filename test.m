clc, clear;
close all

% K = 434;
% taum = 0.67;
% taue = 0.0035;

K = 1.56 * 180 / pi;
taue = 0.0039035;
taum = 0.984871194396488;

bandwidth = 18 * 2 * pi;
wc_max = 650;
[P, G, para] = direct_design(bandwidth, wc_max, K, taum, taue);

figurename('ֱ����ƿ���');
margin(P * G);
grid on

figurename('ֱ����Ʊջ�');
bode(P * G / (1 + P * G));
grid on

Design_Lowgain

%% �Ƿ���Ҫ�����ܷ���Ƴ�������
[trap, later, bfailure, data_check] = wave_repair(P * Glow, G, para, 250, 0, bandwidth);
K = P * G * later.G * Glow;
for i = 1 : trap.num
    K = K * trap.G(i);
end
figurename('�ݲ��˲���');
margin(K);
grid on
figurename('�ݲ��˲����ջ�');
bode(K / (1 + K));
grid on


% Design_Lowgain
% 
% figurename('��Ƶ����');
% K = K * Glow;
% margin(K);
% grid on
% figurename('��Ƶ�ջ�');
% bode(K / (1 + K));
% grid on


autoArrangeFigures;
