clc, clear
close all

%% ����ģ���Ƿ���������
ncount = 20;
data.mag = zeros(ncount, 1);
data.phi = zeros(ncount, 1);
data.fre = linspace(1, 20) * 2 * pi;

for i = 1 : ncount
    data.mag(i) = 1 + 0.02 * (i - 10);
    data.phi(i) = 2 * (i - 10);
end

bandwidth = 20 * 2 * pi;
option.type = 'close-loop';

data_check = CLIndic_check(data, bandwidth, option);
%% ����

G = tf(1, [50, 1]) * 50;
margin(G);
grid on
G_C = G / (1 + G);
figurename('bihuan');
bode(G_C);
grid on
data.mag = zeros(ncount, 1);
data.phi = zeros(ncount, 1);
data.fre = linspace(1, 20) * 2 * pi;
[mag, phi] = bode(G, data.fre);
for i = 1 : ncount
   data.mag(i) = mag(1, 1, i);
   data.phi(i) = phi(1, 1, i);
end
option.type = 'open-loop';
data_check = CLIndic_check(data, bandwidth, option);
%% ���ͨ��
G = tf(100, [1, 5, 100]);
figurename('г����');
bode(G);
grid on
option.type = 'transfer';

[fre, Rp, flag] = Rp_check(G, bandwidth, 5, 5, option);

%% ͨ��
[mag, phi, fre] = bode(G);
data.mag = zeros(ncount, 1);
data.phi = zeros(ncount, 1);
data.fre = fre;
for i = 1 : 1 : length(fre)
   data.mag(i) = mag(1, 1, i);
   data.phi(i) = phi(1, 1, i);
end
option.type = 'discrete';
[fre, Rp, flag] = Rp_check(data, bandwidth, 5, 5, option);
%% ͨ��


option.type = 'transfer';
G = tf(1, [100, 1]) * 50;
figurename('wendingxing');
margin(G);
grid on
[bStable, bGm, bPm, bPhi, bWc] = Stability_check(G, G, 1, 20, 100, 150, 100, 2, option);





autoArrangeFigures;


