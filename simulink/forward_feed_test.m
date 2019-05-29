clc, clear
close all

tau = 0.0005 * 15;

w = 1 : 100;
w = w * 2 * pi;

s = tf('s');
G1 = exp(-tau * s);
[mag1, phi1] = bode(G1, w);
mag1 = 20 * log10(reshape(mag1, [length(w), 1]));
phi1 = reshape(phi1, [length(w), 1]);
% complex_e = complex(cos(-tau * w), sin(-tau * w));
% mag1 = 20 * log10(abs(complex_e));
% phi1 = angle(complex_e) / pi * 180;

G = tf(1, [tau, 1]);
[mag2, phi2] = bode(G, w);
mag2 = 20 * log10(reshape(mag2, [length(w), 1]));
phi2 = reshape(phi2, [length(w), 1]);

figurename('bode');
subplot 211
semilogx(w, mag1, 'r*-');
hold on
grid on
semilogx(w, mag2, 'b*-');
subplot 212
semilogx(w, phi1, 'r*-');
hold on
grid on
semilogx(w, phi2, 'b*-');

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
figurename('a');
bode(G_speed);
grid on
hold on
Gs = tf([1 0],1);
bode(Gs);

