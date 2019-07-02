clc, clear
close all

load('����ESO10hz�ջ�.csv');
load('����ESO15hz�ջ�.csv');
X15QC = load('����ESO15hzǰ��.csv');
X15ZL = load('����ESO15hzָ��Ԥ����.csv');


fre = linspace(1, 10, 10) * 2 * pi;
mag1 = 20 * log10(X__ESO10hz__(:, 1));
phi1 = X__ESO10hz__(:, 2);
mag2 = 20 * log10(X__ESO10hz__(:, 3));
phi2 = X__ESO10hz__(:, 4);
mag3 = 20 * log10(X__ESO10hz__(:, 5));
phi3 = X__ESO10hz__(:, 6);
mag4 = 20 * log10(X__ESO10hz__(:, 7));
phi4 = X__ESO10hz__(:, 8);


figurename('10hz');
subplot 211
semilogx(fre, mag1, 'bo-');
hold on
grid on
% semilogx(fre, mag2, 'ro-');
semilogx(fre, mag3, 'g*-');
semilogx(fre, mag4, 'r.-');
% legend('����', '����', '����', '����');
legend('����', '����', '����');

ylabel('��ֵ(db)');

subplot 212
semilogx(fre, phi1, 'bo-');
hold on
grid on
% semilogx(fre, phi2, 'ro-');
semilogx(fre, phi3, 'g*-');
semilogx(fre, phi4, 'r.-');
ylabel('��λ(��)');
xlabel('��Ƶ��(rad/s)');

figurename('10hz�Ա�');
subplot 211
% semilogx(fre, mag1, 'bo-');
hold on
grid on
semilogx(fre, mag3 - mag1, 'g*-');
semilogx(fre, mag4 - mag1, 'b.-');
% legend('����', '����', '����', '����');
legend('����', '����');
ylabel('��ֵ��(db)');

subplot 212
% semilogx(fre, phi1, 'bo-');
hold on
grid on
semilogx(fre, phi3 - phi1, 'g*-');
semilogx(fre, phi4 - phi1, 'b.-');
ylabel('��λ��(��)');
xlabel('��Ƶ��(rad/s)');


fre = linspace(1, 15, 15) * 2 * pi;
mag1 = 20 * log10(X__ESO15hz__(:, 1));
phi1 = X__ESO15hz__(:, 2);
mag2 = 20 * log10(X__ESO15hz__(:, 3));
phi2 = X__ESO15hz__(:, 4);
mag3 = 20 * log10(X__ESO15hz__(:, 5));
phi3 = X__ESO15hz__(:, 6);
mag4 = 20 * log10(X__ESO15hz__(:, 7));
phi4 = X__ESO15hz__(:, 8);

mag11 = mag1;
phi11 = phi1;

figurename('15hz');
subplot 211
semilogx(fre, mag1, 'bo-');
hold on
grid on
% semilogx(fre, mag2, 'ro-');
semilogx(fre, mag3, 'g*-');
semilogx(fre, mag4, 'r.-');
% legend('����', '����', '����', '����');
legend('����', '����', '����');

ylabel('��ֵ(db)');

subplot 212
semilogx(fre, phi1, 'bo-');
hold on
grid on
% semilogx(fre, phi2, 'ro-');
semilogx(fre, phi3, 'g*-');
semilogx(fre, phi4, 'r.-');
ylabel('��λ(��)');
xlabel('��Ƶ��(rad/s)');

figurename('15hz�Ա�');
subplot 211
% semilogx(fre, mag1, 'bo-');
hold on
grid on
semilogx(fre, mag3 - mag1, 'g*-');
semilogx(fre, mag4 - mag1, 'b.-');
% legend('����', '����', '����', '����');
legend('����', '����');
ylabel('��ֵ��(db)');

subplot 212
% semilogx(fre, phi1, 'bo-');
hold on
grid on
semilogx(fre, phi3 - phi1, 'g*-');
semilogx(fre, phi4 - phi1, 'b.-');
ylabel('��λ��(��)');
xlabel('��Ƶ��(rad/s)');

fre = linspace(1, 15, 15) * 2 * pi;
mag1 = 20 * log10(X15QC(:, 1));
phi1 = X15QC(:, 2);
mag2 = 20 * log10(X15QC(:, 3));
phi2 = X15QC(:, 4);
mag3 = 20 * log10(X15QC(:, 5));
phi3 = X15QC(:, 6);
mag4 = 20 * log10(X15QC(:, 7));
phi4 = X15QC(:, 8);

mag111 = mag1;
phi111 = phi1;

figurename('15hzǰ��');
subplot 211
semilogx(fre, mag1, 'b*-');
hold on
grid on
semilogx(fre, mag2, 'ro-');
semilogx(fre, mag11, 'g.-');


% semilogx(fre, mag3, 'g*-');
% semilogx(fre, mag4, 'r.-');
% legend('����', '����', '����', '����');
% legend('����', '����', '����');
legend('ʹ��ǰ��', 'ǰ������', 'δʹ��ǰ��');


ylabel('��ֵ(db)');

subplot 212
semilogx(fre, phi1, 'b*-');
hold on
grid on
semilogx(fre, phi2, 'ro-');
semilogx(fre, phi11, 'g.-');

% semilogx(fre, phi3, 'g*-');
% semilogx(fre, phi4, 'r.-');
ylabel('��λ(��)');
xlabel('��Ƶ��(rad/s)');

fre = linspace(1, 15, 15) * 2 * pi;
mag1 = 20 * log10(X15ZL(:, 1));
phi1 = X15ZL(:, 2);
mag2 = X15ZL(:, 3);
phi2 = X15ZL(:, 4);
mag3 = 20 * log10(X15ZL(:, 5));
phi3 = X15ZL(:, 6);
mag4 = 20 * log10(X15ZL(:, 7));
phi4 = X15ZL(:, 8);

figurename('15hzԤ����');
subplot 211
semilogx(fre, mag1, 'b*-');
hold on
grid on
semilogx(fre, mag2, 'ro-');
semilogx(fre, mag111, 'g.-');

% semilogx(fre, mag3, 'g*-');
% semilogx(fre, mag4, 'r.-');
% legend('����', '����', '����', '����');legend('����', '����', '����');
% legend('����', '����', '����');
legend('ʹ��ָ��Ԥ����', 'ָ��Ԥ��������', 'δʹ��ָ��Ԥ����');

ylabel('��ֵ(db)');

subplot 212
semilogx(fre, phi1, 'b*-');
hold on
grid on
semilogx(fre, phi2, 'ro-');
semilogx(fre, phi111, 'g.-');

% semilogx(fre, phi3, 'g*-');
% semilogx(fre, phi4, 'r.-');
ylabel('��λ(��)');
xlabel('��Ƶ��(rad/s)');

load('opensweepcompare.csv');
mag1 = 20 * log10(opensweepcompare(:, 1));
phi1 = opensweepcompare(:, 2);
mag2 = 20 * log10(opensweepcompare(:, 3));
phi2 = opensweepcompare(:, 4);
mag3 = 20 * log10(opensweepcompare(:, 5));
phi3 = opensweepcompare(:, 6);

figurename('esoopen');
subplot 211
fre = linspace(1, 50, 50) * 2 * pi;
semilogx(fre, mag1, 'r.-');
hold on
semilogx(fre, mag2, 'bo-');
semilogx(fre, mag3, 'g*-');
grid on
legend('����', '����', '����');
ylabel('��ֵ(db)');

subplot 212 
semilogx(fre, phi1, 'r.-');
hold on
semilogx(fre, phi2, 'bo-');
semilogx(fre, phi3, 'g*-');
grid on
ylabel('��λ(��)');
xlabel('��Ƶ��(rad/s)');

figurename('����');
subplot 211
fre = linspace(1, 50, 50) * 2 * pi;
semilogx(fre, mag1 - mag2, 'b.-');
hold on
% semilogx(fre, mag2, 'bo-');
semilogx(fre, mag3 - mag2, 'g*-');
grid on
% legend('����', '����', '����');
legend('����', '����');

ylabel('��ֵ��(db)');

subplot 212 
semilogx(fre, phi1 - phi2, 'b.-');
hold on
% semilogx(fre, phi2, 'bo-');
semilogx(fre, phi3 - phi2, 'g*-');
grid on
ylabel('��λ��(��)');
xlabel('��Ƶ��(rad/s)');



