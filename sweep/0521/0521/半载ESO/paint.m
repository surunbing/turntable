clc, clear
close all

load('10hz.csv');
load('15hz.csv');

fre = linspace(1, 10, 10) * 2 * pi;
mag1 = 20 * log10(X10hz(:, 1));
phi1 = X10hz(:, 2);
mag2 = 20 * log10(X10hz(:, 3));
phi2 = X10hz(:, 4);


figurename('10hz');
subplot 211
semilogx(fre, mag1, 'b*-');
hold on
grid on
semilogx(fre, mag2, 'ro-');
legend('实际', '理论');
ylabel('幅值(db)');

subplot 212
semilogx(fre, phi1, 'b*-');
hold on
grid on
semilogx(fre, phi2, 'ro-');
ylabel('相位(°)');
xlabel('角频率(rad/s)');


fre = linspace(1, 15, 15) * 2 * pi;
mag1 = 20 * log10(X15hz(:, 1));
phi1 = X15hz(:, 2);
mag2 = 20 * log10(X15hz(:, 3));
phi2 = X15hz(:, 4);


figurename('15hz');
subplot 211
semilogx(fre, mag1, 'b*-');
hold on
grid on
semilogx(fre, mag2, 'ro-');
legend('实际', '理论');
ylabel('幅值(db)');

subplot 212
semilogx(fre, phi1, 'b*-');
hold on
grid on
semilogx(fre, phi2, 'ro-');
ylabel('相位(°)');
xlabel('角频率(rad/s)');
