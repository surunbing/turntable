clc, clear
close all
e = 0.92;
T = 4;
f = 16 * 2 * pi;
e1 = 0.92;
T1 = 100;
f1 = 10 * 2 * pi;
trap = tf([1 e * T f * f], [1, T f * f]);
trap1 = tf([1 e1 * T1 f1 * f1], [1, T1 f1 * f1]);
bode(trap);
hold on
bode(trap1);
hold on
bode(trap * trap1);
grid on