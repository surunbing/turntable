clc, clear
close all

load('半载ESO10hz闭环.csv');
load('半载ESO15hz闭环.csv');
X15QC = load('半载ESO15hz前馈.csv');
X15ZL = load('半载ESO15hz指令预处理.csv');


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
% legend('半载', '理论', '满载', '空载');
legend('半载', '满载', '空载');

ylabel('幅值(db)');

subplot 212
semilogx(fre, phi1, 'bo-');
hold on
grid on
% semilogx(fre, phi2, 'ro-');
semilogx(fre, phi3, 'g*-');
semilogx(fre, phi4, 'r.-');
ylabel('相位(°)');
xlabel('角频率(rad/s)');

figurename('10hz对比');
subplot 211
% semilogx(fre, mag1, 'bo-');
hold on
grid on
semilogx(fre, mag3 - mag1, 'g*-');
semilogx(fre, mag4 - mag1, 'b.-');
% legend('半载', '理论', '满载', '空载');
legend('满载', '空载');
ylabel('幅值差(db)');

subplot 212
% semilogx(fre, phi1, 'bo-');
hold on
grid on
semilogx(fre, phi3 - phi1, 'g*-');
semilogx(fre, phi4 - phi1, 'b.-');
ylabel('相位差(°)');
xlabel('角频率(rad/s)');


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
% legend('半载', '理论', '满载', '空载');
legend('半载', '满载', '空载');

ylabel('幅值(db)');

subplot 212
semilogx(fre, phi1, 'bo-');
hold on
grid on
% semilogx(fre, phi2, 'ro-');
semilogx(fre, phi3, 'g*-');
semilogx(fre, phi4, 'r.-');
ylabel('相位(°)');
xlabel('角频率(rad/s)');

figurename('15hz对比');
subplot 211
% semilogx(fre, mag1, 'bo-');
hold on
grid on
semilogx(fre, mag3 - mag1, 'g*-');
semilogx(fre, mag4 - mag1, 'b.-');
% legend('半载', '理论', '满载', '空载');
legend('满载', '空载');
ylabel('幅值差(db)');

subplot 212
% semilogx(fre, phi1, 'bo-');
hold on
grid on
semilogx(fre, phi3 - phi1, 'g*-');
semilogx(fre, phi4 - phi1, 'b.-');
ylabel('相位差(°)');
xlabel('角频率(rad/s)');

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

figurename('15hz前馈');
subplot 211
semilogx(fre, mag1, 'b*-');
hold on
grid on
semilogx(fre, mag2, 'ro-');
semilogx(fre, mag11, 'g.-');


% semilogx(fre, mag3, 'g*-');
% semilogx(fre, mag4, 'r.-');
% legend('半载', '理论', '满载', '空载');
% legend('半载', '满载', '空载');
legend('使用前馈', '前馈理论', '未使用前馈');


ylabel('幅值(db)');

subplot 212
semilogx(fre, phi1, 'b*-');
hold on
grid on
semilogx(fre, phi2, 'ro-');
semilogx(fre, phi11, 'g.-');

% semilogx(fre, phi3, 'g*-');
% semilogx(fre, phi4, 'r.-');
ylabel('相位(°)');
xlabel('角频率(rad/s)');

fre = linspace(1, 15, 15) * 2 * pi;
mag1 = 20 * log10(X15ZL(:, 1));
phi1 = X15ZL(:, 2);
mag2 = X15ZL(:, 3);
phi2 = X15ZL(:, 4);
mag3 = 20 * log10(X15ZL(:, 5));
phi3 = X15ZL(:, 6);
mag4 = 20 * log10(X15ZL(:, 7));
phi4 = X15ZL(:, 8);

figurename('15hz预处理');
subplot 211
semilogx(fre, mag1, 'b*-');
hold on
grid on
semilogx(fre, mag2, 'ro-');
semilogx(fre, mag111, 'g.-');

% semilogx(fre, mag3, 'g*-');
% semilogx(fre, mag4, 'r.-');
% legend('半载', '理论', '满载', '空载');legend('半载', '满载', '空载');
% legend('半载', '满载', '空载');
legend('使用指令预处理', '指令预处理理论', '未使用指令预处理');

ylabel('幅值(db)');

subplot 212
semilogx(fre, phi1, 'b*-');
hold on
grid on
semilogx(fre, phi2, 'ro-');
semilogx(fre, phi111, 'g.-');

% semilogx(fre, phi3, 'g*-');
% semilogx(fre, phi4, 'r.-');
ylabel('相位(°)');
xlabel('角频率(rad/s)');

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
legend('空载', '半载', '满载');
ylabel('幅值(db)');

subplot 212 
semilogx(fre, phi1, 'r.-');
hold on
semilogx(fre, phi2, 'bo-');
semilogx(fre, phi3, 'g*-');
grid on
ylabel('相位(°)');
xlabel('角频率(rad/s)');

figurename('差异');
subplot 211
fre = linspace(1, 50, 50) * 2 * pi;
semilogx(fre, mag1 - mag2, 'b.-');
hold on
% semilogx(fre, mag2, 'bo-');
semilogx(fre, mag3 - mag2, 'g*-');
grid on
% legend('空载', '半载', '满载');
legend('空载', '满载');

ylabel('幅值差(db)');

subplot 212 
semilogx(fre, phi1 - phi2, 'b.-');
hold on
% semilogx(fre, phi2, 'bo-');
semilogx(fre, phi3 - phi2, 'g*-');
grid on
ylabel('相位差(°)');
xlabel('角频率(rad/s)');



