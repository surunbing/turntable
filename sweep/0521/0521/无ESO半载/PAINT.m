clc, clear
close all

load('10hz.csv');
load('15hz.csv');
QK = load('15hz ǰ��.csv');
ZL = load('15hzָ��Ԥ����.csv');
DIR = load('dir.csv');

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
legend('ʵ��', '����');
ylabel('��ֵ(db)');

subplot 212
semilogx(fre, phi1, 'b*-');
hold on
grid on
semilogx(fre, phi2, 'ro-');
ylabel('��λ(��)');
xlabel('��Ƶ��(rad/s)');


fre = linspace(1, 15, 15) * 2 * pi;
mag1 = 20 * log10(X15hz(:, 1));
phi1 = X15hz(:, 2);
mag2 = 20 * log10(X15hz(:, 3));
phi2 = X15hz(:, 4);

mag11 = mag1;
phi11 = phi1;


figurename('15hz');
subplot 211
semilogx(fre, mag1, 'b*-');
hold on
grid on
semilogx(fre, mag2, 'ro-');
legend('ʵ��', '����');
ylabel('��ֵ(db)');

subplot 212
semilogx(fre, phi1, 'b*-');
hold on
grid on
semilogx(fre, phi2, 'ro-');
ylabel('��λ(��)');

magdiff = mag1 - mag2;
phidiff = phi1 - phi2;

fre = linspace(1, 15, 15) * 2 * pi;
mag1 = 20 * log10(QK(:, 1)) + magdiff;
phi1 = QK(:, 2) + phidiff;
mag2 = 20 * log10(QK(:, 3)) ;
phi2 = QK(:, 4);

mag111 = mag1;
phi111 = phi1;


figurename('15hzǰ��');
subplot 211
semilogx(fre, mag1, 'b*-');
hold on
grid on
semilogx(fre, mag2, 'ro-');
semilogx(fre, mag11, 'g.-');
legend('ʹ��ǰ��', 'ǰ������', 'δʹ��ǰ��');
ylabel('��ֵ(db)');

subplot 212
semilogx(fre, phi1, 'b*-');
hold on
grid on
semilogx(fre, phi2, 'ro-');
semilogx(fre, phi11, 'g.-');
ylabel('��λ(��)');
xlabel('��Ƶ��(rad/s)');
xlabel('��Ƶ��(rad/s)');

fre = linspace(1, 15, 15) * 2 * pi;
mag1 = 20 * log10(ZL(:, 1)) + magdiff;
phi1 = ZL(:, 2) + phidiff;
mag2 = 20 * log10(ZL(:, 3));
phi2 = ZL(:, 4);


figurename('15hzָ��');
subplot 211
semilogx(fre, mag1, 'b*-');
hold on
grid on
semilogx(fre, mag2, 'ro-');
semilogx(fre, mag111, 'g.-');
legend('ʹ��ָ��Ԥ����', 'ָ��Ԥ��������', 'δʹ��ָ��Ԥ����');
ylabel('��ֵ(db)');

subplot 212
semilogx(fre, phi1, 'b*-');
hold on
grid on
semilogx(fre, phi2, 'ro-');
semilogx(fre, phi111, 'g.-');
ylabel('��λ(��)');
xlabel('��Ƶ��(rad/s)');
xlabel('��Ƶ��(rad/s)');


fre = linspace(1, 10, 10) * 2 * pi;
mag1 = 20 * log10(DIR(:, 1));
phi1 = DIR(:, 2);
mag2 = 20 * log10(DIR(:, 3));
phi2 = DIR(:, 4);


figurename('DIR');
subplot 211
semilogx(fre, mag1, 'b*-');
hold on
grid on
semilogx(fre, mag2, 'ro-');
legend('ʵ��', '����');
ylabel('��ֵ(db)');

subplot 212
semilogx(fre, phi1, 'b*-');
hold on
grid on
semilogx(fre, phi2, 'ro-');
ylabel('��λ(��)');
xlabel('��Ƶ��(rad/s)');

