clc, clear
close all

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

%% 载入相关数据对比
%% 载入开环数据对比
data_eso = load('0516sweepeso-3.csv');
data_noeso = load('0516sweepnoeso.csv');
eso_open.fre = data_eso(:, 1) * 2 * pi;
eso_open.mag = 20 * log10(data_eso(:, 2));
eso_open.phi = data_eso(:, 3);
eso_open.complex = data_eso(:, 2) .* complex(cos(data_eso(:, 3) / 180 * pi), sin(data_eso(:, 3) / 180 * pi));
noeso_open.fre = data_noeso(:, 1) * 2 * pi;
noeso_open.mag = 20 * log10(data_noeso(:, 2));
noeso_open.phi = data_noeso(:, 3);
noeso_open.complex = data_noeso(:, 2) .* complex(cos(data_noeso(:, 3) / 180 * pi), sin(data_noeso(:, 3) / 180 * pi));

figurename('eso扫频对比');
subplot 211
semilogx(eso_open.fre, eso_open.mag, 'r*-');
hold on
grid on
semilogx(noeso_open.fre, noeso_open.mag, 'b*-');

subplot 212
semilogx(eso_open.fre, eso_open.phi, 'r*-');
hold on
grid on
semilogx(noeso_open.fre, noeso_open.phi, 'b*-');

legend('eso', '原始');


%% 载入数据后，对比理想的闭环特性
%% 计算控制器
T = 1 / (10 * 2 * pi) / 10; 
Inertial = tf(1, [T, 1]);

%% 30穿越  20 * 3
fre = 30 * 2 * pi;
aangle = 20;
m = sin((aangle) / 180 * pi);
alpha = (1 - m) / (1 + m);
tau = 1 / (sqrt(alpha) * fre);
T = alpha * tau;
G_advance1 = tf([tau, 1], [T, 1]);

%% 30穿越  10 * 1
fre = 30 * 2 * pi;
aangle = 10;
m = sin((aangle) / 180 * pi);
alpha = (1 - m) / (1 + m);
tau = 1 / (sqrt(alpha) * fre);
T = alpha * tau;
G_advance2 = tf([tau, 1], [T, 1]);

%% 30穿越  15 * 1
fre = 30 * 2 * pi;
aangle = 15;
m = sin((aangle) / 180 * pi);
alpha = (1 - m) / (1 + m);
tau = 1 / (sqrt(alpha) * fre);
T = alpha * tau;
G_advance3 = tf([tau, 1], [T, 1]);

%% 可以对比迟后
alpha = 3;
fre = 1 * 2 * pi;
tau = 1 / (sqrt(alpha * fre));
G_later1 = tf([tau, 1], [alpha * tau, 1]);


gain1 = 1.1;
gain2 = 51.5464;

controller_noeso.G = Inertial * G_advance1 * G_advance1 * G_advance1 * G_advance2 * G_advance3 * gain2;
controller_eso.G = Inertial * G_advance1 * G_advance1 * G_advance1 * G_advance2 * G_advance3 * gain1 * gain2;

[mag, phi] = bode(controller_noeso.G, noeso_open.fre);
controller_noeso.mag = 20 * log10(reshape(mag, [length(noeso_open.fre), 1]));
controller_noeso.phi = reshape(phi, [length(noeso_open.fre), 1]);
controller_noeso.complex = reshape(mag, [length(noeso_open.fre), 1]) .* complex(cos(controller_noeso.phi / 180 * pi), sin(controller_noeso.phi / 180 * pi));

[mag, phi] = bode(controller_eso.G, eso_open.fre);
controller_eso.mag = 20 * log10(reshape(mag, [length(eso_open.fre), 1]));
controller_eso.phi = reshape(phi, [length(eso_open.fre), 1]);
controller_eso.complex = reshape(mag, [length(eso_open.fre), 1]) .* complex(cos(controller_eso.phi / 180 * pi), sin(controller_eso.phi / 180 * pi));

open_forward_noeso = noeso_open.complex .* controller_noeso.complex;
open_forward_eso = eso_open.complex .* controller_eso.complex;

figurename('前向通道对比');
subplot 211
semilogx(noeso_open.fre, 20 * log10(abs(open_forward_noeso)), 'r*-');
grid on
hold on
semilogx(eso_open.fre, 20 * log10(abs(open_forward_eso)), 'b*-');

subplot 212
semilogx(noeso_open.fre, angle(open_forward_noeso) / pi * 180, 'r*-');
grid on
hold on
semilogx(eso_open.fre, angle(open_forward_eso) / pi * 180, 'b*-');

legend('noeso', 'eso');

%% 理想闭环与实际闭环
close_forward_noeso = open_forward_noeso ./ (1 + open_forward_noeso);
close_forward_eso = open_forward_eso ./ (1 + open_forward_eso);

figurename('闭环通道理论');
subplot 211
semilogx(noeso_open.fre, 20 * log10(abs(close_forward_noeso)), 'r*-');
grid on
hold on
semilogx(eso_open.fre, 20 * log10(abs(close_forward_eso)), 'b*-');

subplot 212
semilogx(noeso_open.fre, angle(close_forward_noeso) / pi * 180, 'r*-');
grid on
hold on
semilogx(eso_open.fre, angle(close_forward_eso) / pi * 180, 'b*-');

legend('noeso', 'eso');

figurename('闭环通道实际');
data_noeso_close = load('0516closenoeso.csv');
data_eso_close = load('0516closeeso-gain11.csv');

noeso_close.fre = data_noeso_close(:, 1) * 2 * pi;
noeso_close.mag = 20 * log10(data_noeso_close(:, 2));
noeso_close.phi = data_noeso_close(:, 3);

eso_close.fre = data_eso_close(:, 1) * 2 * pi;
eso_close.mag = 20 * log10(data_eso_close(:, 2));
eso_close.phi = data_eso_close(:, 3);

subplot 211
semilogx(noeso_close.fre, noeso_close.mag, 'r*-');
grid on
hold on
semilogx(eso_close.fre, eso_close.mag, 'b*-');

subplot 212
semilogx(noeso_close.fre, noeso_close.phi, 'r*-');
grid on
hold on
semilogx(eso_close.fre, eso_close.phi, 'b*-');

legend('noeso', 'eso');

%% 对比 没有eso
figurename('无eso对比');
subplot 211
semilogx(noeso_open.fre, 20 * log10(abs(close_forward_noeso)), 'r*-');
grid on
hold on
semilogx(noeso_close.fre, noeso_close.mag, 'b*-');

subplot 212
semilogx(noeso_open.fre,  angle(close_forward_noeso) / pi * 180, 'r*-');
grid on
hold on
semilogx(noeso_close.fre, noeso_close.phi, 'b*-');

legend('理论', '实际');

%% 对比 有eso
figurename('eso对比');
subplot 211
semilogx(eso_open.fre, 20 * log10(abs(close_forward_eso)), 'r*-');
grid on
hold on
semilogx(eso_close.fre, eso_close.mag, 'b*-');

subplot 212
semilogx(eso_open.fre,  angle(close_forward_eso) / pi * 180, 'r*-');
grid on
hold on
semilogx(eso_close.fre, eso_close.phi, 'b*-');

legend('理论', '实际');


%% 前馈验证
K = 355.1054;
taue = 0.001668255197954;
taum = 2.196498493001917;
% G_model = tf(K, [taue * taum taue + taum 1 0]);
s = tf('s');
G_auxiliary = tf(1, conv([1 / (200 * 2 * pi), 1], [1 / (200 * 2 * pi), 1]));
forward_app = tf(K, [taue * taum taue + taum 1]);
forward.G = 0.5 / forward_app * G_auxiliary * G_speed;
[mag, phi] = bode(forward.G, eso_open.fre);
forward.mag = reshape(mag, [length(eso_open.fre), 1]);
forward.phi = reshape(phi, [length(eso_open.fre), 1]);
forward.complex = forward.mag .* complex(cos(forward.phi / 180 * pi), sin(forward.phi / 180 * pi));

% 理论noeso
close_noeso_forward = (open_forward_noeso + noeso_open.complex .* forward.complex) ./ (1 + open_forward_noeso);
close_eso_forward = (open_forward_eso + eso_open.complex .* forward.complex) ./ (1 + open_forward_eso);

data_close_noeso_forward = load('noesosweepclose-for0.5-0516.csv');
data_close_eso_forward = load('esosweepclose-for0.5-0516.csv');

figurename('noeso对比前馈');
subplot 211
semilogx(eso_open.fre, 20 * log10(abs(close_noeso_forward)), 'r*-');
grid on
hold on
semilogx(data_close_noeso_forward(:, 1) * 2 * pi, 20 * log10(data_close_noeso_forward(:, 2)), 'b*-');

subplot 212
semilogx(eso_open.fre,  angle(close_noeso_forward) / pi * 180, 'r*-');
grid on
hold on
semilogx(data_close_noeso_forward(:, 1) * 2 * pi, data_close_noeso_forward(:, 3), 'b*-');

legend('理论', '实际');

figurename('eso对比前馈');
subplot 211
semilogx(eso_open.fre, 20 * log10(abs(close_eso_forward)), 'r*-');
grid on
hold on
semilogx(data_close_eso_forward(:, 1) * 2 * pi, 20 * log10(data_close_eso_forward(:, 2)), 'b*-');

subplot 212
semilogx(eso_open.fre,  angle(close_eso_forward) / pi * 180, 'r*-');
grid on
hold on
semilogx(data_close_eso_forward(:, 1) * 2 * pi, data_close_eso_forward(:, 3), 'b*-');

legend('理论', '实际');

%% 迟后验证
alpha = 3;
fre = 1 * 2 * pi;
tau = 1 / (sqrt(alpha * fre));
G_later3 = alpha * tf([tau, 1], [alpha * tau, 1]);

gain3 = 3;

[mag, phi] = bode(G_later3, eso_open.fre);
later.mag = reshape(mag, [length(eso_open.fre), 1]);
later.phi = reshape(phi, [length(eso_open.fre), 1]);
later.complex = later.mag .* complex(cos(later.phi / 180 * pi), sin(later.phi / 180 * pi));

open_forward_noeso_later = noeso_open.complex .* controller_noeso.complex .* later.complex;
open_forward_eso_later = eso_open.complex .* controller_eso.complex .* later.complex;

close_forward_noeso_later = open_forward_noeso_later ./ (1 + open_forward_noeso_later);
close_forward_eso_later = open_forward_eso_later ./ (1 + open_forward_eso_later);

figurename('闭环通道理论later-eso');
subplot 211
semilogx(eso_open.fre, 20 * log10(abs(close_forward_noeso)), 'r*-');
grid on
hold on
semilogx(eso_open.fre, 20 * log10(abs(close_forward_noeso_later)), 'b*-');

subplot 212
semilogx(eso_open.fre, angle(close_forward_noeso) / pi * 180, 'r*-');
grid on
hold on
semilogx(eso_open.fre, angle(close_forward_noeso_later) / pi * 180, 'b*-');

legend('nolater', 'later');

autoArrangeFigures;


% K = 355.1054;
% 
% taue = 0.001668255197954;
% taum = 2.196498493001917;
% 
% G_model = tf(K, [taue * taum taue + taum 1 0]);
% 
% T = 1 / (10 * 2 * pi) / 10; 
% Inertial = tf(1, [T, 1]);
% 
% alpha = 3;
% fre = 0.1 * 2 * pi;
% tau = 1 / (sqrt(alpha * fre));
% G_later1 = tf([tau, 1], [alpha * tau, 1]);
% 
% alpha = 3;
% fre = 0.5 * 2 * pi;
% tau = 1 / (sqrt(alpha * fre));
% G_later2 = tf([tau, 1], [alpha * tau, 1]);
% 
% alpha = 3;
% fre = 1 * 2 * pi;
% tau = 1 / (sqrt(alpha * fre));
% G_later3 = tf([tau, 1], [alpha * tau, 1]);
% 
% frequence = linspace(1, 15, 15) * 2 * pi;
% [mag, phi] = bode((P * G + G * forward)/ (1 + P * G), frequence);
% % [mag, phi] = bode((P * G)/(1 + P * G), frequence);
% mag = reshape(mag, [15, 1]);
% phi = reshape(phi, [15, 1]);
% phi = phi;



% a = [1, 0];
% b = [0.001, 1];
% TSp = 0.0005;
% [dNumf, dDenf] = c2dm(a, b, TSp, 'tustin');
% 
% ts = 0.0005;
% G = 1 / ts / 100 * (1 + tf(1, [ts, 1])  + tf(1, [2 * ts, 1])  + tf(1, [3 * ts, 1]) + tf(1, [4 * ts, 1]) + tf(1, [5 * ts, 1]) + tf(1, [6 * ts, 1]) + tf(1, [7 * ts, 1]) + tf(1, [8 * ts, 1]) + tf(1, [9 * ts, 1])- tf(1, [10 * ts, 1]) - tf(1, [11 * ts, 1]) -tf(1, [12 * ts, 1]) - tf(1, [13 * ts, 1]) - tf(1, [14 * ts, 1]) - tf(1, [15 * ts, 1]) - tf(1, [16 * ts, 1]) - tf(1, [17 * ts, 1]) - tf(1, [18 * ts, 1]) - tf(1, [19 * ts, 1]));
% 
% bode(G);
% grid on
    
