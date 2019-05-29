clc, clear
close all

K = 1.56 * 180 / pi;
taue = 0.0039035;
taum = 0.984871194396488;

taum2 = taum * 1;

G =  tf(K, [taue * taum2, taue + taum2, 1, 0]);
G2 = tf([taue * taum, taue + taum, 1, 0], K);


wq = 60 * 2 * pi;
% beta4 = 0.2756;
% beta1 = 0.9528;
% beta2 = 1.4539;
% beta3 = 0.7426;

beta4 = 1;
beta1 = 4;
beta2 = 6;
beta3 = 4;
Q = tf(beta4 * (wq ^ 4), [1, beta1 * wq, beta2 * wq * wq, beta3 * wq * wq * wq, beta4 * wq * wq * wq * wq]);
figurename('ESO_filter');
bode(Q);
grid on

G_eso = (G2 * G - 1) * Q;
GESO_CLOSE = 1 / (1 + G_eso);

load('controller_10_3_nolow.mat');

figurename('close');
G_close = (P * GESO_CLOSE * G) / (1 + P * GESO_CLOSE * G);
bode(G_close);
grid on


ts = 0.0005;
forward_n = 15;
G_speed = 1 / ts / forward_n / forward_n;
G_speed_sum = 0;
s = tf('s');
for i = 1 : forward_n
    G_speed_sum = G_speed_sum + exp(-(i - 1) * ts * s) - exp(-(i + forward_n - 1) * ts * s);%tf(1, [(i - 1) * ts, 1]) - tf(1, [(i + 9) * ts, 1]);
end
% G_speed = 1 / ts / 100 * (1 + tf(1, [ts, 1])  + tf(1, [2 * ts, 1])  + tf(1, [3 * ts, 1]) + tf(1, [4 * ts, 1]) + tf(1, [5 * ts, 1]) + tf(1, [6 * ts, 1]) + tf(1, [7 * ts, 1]) + tf(1, [8 * ts, 1]) + tf(1, [9 * ts, 1])- tf(1, [10 * ts, 1]) - tf(1, [11 * ts, 1]) -tf(1, [12 * ts, 1]) - tf(1, [13 * ts, 1]) - tf(1, [14 * ts, 1]) - tf(1, [15 * ts, 1]) - tf(1, [16 * ts, 1]) - tf(1, [17 * ts, 1]) - tf(1, [18 * ts, 1]) - tf(1, [19 * ts, 1]));
G_speed = G_speed * G_speed_sum;

G_for = tf([taue * taum, taue + taum, 1], K);
G_close = G * P / (1 + G * P);
G_forward = (G * P + 0.5) / (1 + G * P);
G_forward2 = G_for * G_speed * G;
G_for2 = (G * P + 0.5 * G_forward2) / (1 + G * P);

[mag, phi] = bode(G_for * G_speed * G * 0.5, 10 * 2 * pi);
com = mag * complex(cos(phi / 180 * pi), sin(phi / 180 * pi));
[mag1, phi1] = bode( G * P, 10 * 2 * pi);
com1 = mag1 * complex(cos(phi1 / 180 * pi), sin(phi1 / 180 * pi));

angle1 = angle((com + com1) / (1 + com1)) / pi * 180;
angle2 = angle((com1) / (1 + com1)) / pi * 180;
angle3 = angle((com1 + 0.5) / (1 + com1)) / pi * 180;

angle4 = angle((com + com1)) / pi * 180;
angle5 = angle((com1 + 0.5)) / pi * 180;
angle6 = angle(com1) / pi * 180;

figurename('qiankui');
bode(G_forward);
hold on
grid on
bode(G_close);
bode(G_for2);
legend('Ç°À¡','yuanshi', 'chafen');


autoArrangeFigures;
