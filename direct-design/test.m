clc, clear
close all

%% ≤‚ ‘Ω≈±æ
omegan = 405.4322;
xi = 0.8; %0.2807;
x = [xi, omegan];
T = 0.0014;

kg = 20 * log10(GetGm(x, T));
wc = GetWc(x, T);
pm = GetPm(x, T);
Mr = 20 * log10(GetRp(x, T));
dt = GetDt(x, T);

G = tf(omegan * omegan, conv([1, 2 * xi * omegan, omegan * omegan], [T, 1]));
bode(G);
grid on

%% ≤‚ ‘Õ®π˝