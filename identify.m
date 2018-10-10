clc, clear
close all

K = 1.56 * 180 / pi;
taue = 0.0039035;
taum = 0.984871194396488;
G = tf(K, [taue * taum, taum, 1, 0]);
% T = 1 / (bandwidth / 2 / pi) / 10; 
Inertial = 1;%tf(1, [T, 1]);

G = G * Inertial;

figurename('×ªÌ¨²¨ÌØÍ¼');
bode(G);
grid on

time = 20;
nCount = 20 * 2000;
tseries = zeros(nCount, 1);
Input = zeros(nCount, 1);

for i = 1 : nCount
    tseries(i) = i * 0.0005;
    Input(i) = 1;
end

Output = lsim(G, Input, tseries);

figurename('Í¼');
plot(tseries, Input, 'r-');
grid on
hold on
plot(tseries, Output, 'b-');

autoArrangeFigures

