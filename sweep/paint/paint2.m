clc, clear
close all

load('esosweepbanzai.csv');
load('middle_noeso.csv');




K = 3.491260430644539e+02;     %% ET205半载
taue = 0.001726113287549;
taum = 2.465408813244642;

ncount = 50;
G_model = tf(K, [taue * taum * 1 taue + taum * 1 1 0]);
data.fre = linspace(1, ncount, ncount)' * 2 * pi;
[mag, phi] = bode(G_model, data.fre);
data.mag = 20 .* log10(reshape(mag, [ncount, 1]));
data.phi = reshape(phi, [ncount, 1]);

mag1 = 20 * log10(esosweepbanzai(:, 2));
phi1 = esosweepbanzai(:, 3);

mag2 = 20 * log10(middle_noeso(:, 2));
phi2 = middle_noeso(:, 3);

figurename('??');
subplot 211
semilogx(data.fre, data.mag, 'ro-');
hold on
grid on
semilogx(data.fre, mag1, 'b*-');
semilogx(data.fre, mag2, 'g*-');
ylabel('幅值(db)');
legend('理论', 'ESO', '无ES0');

subplot 212
semilogx(data.fre, data.phi, 'ro-');
hold on
grid on
semilogx(data.fre, phi1, 'b*-');
semilogx(data.fre, phi2, 'g*-');
ylabel('幅值(db)');
xlabel('角频率(rad/s)');

figurename('??');
subplot 211
hold on
grid on
semilogx(data.fre,data.mag- mag1, 'b*-');
semilogx(data.fre, data.mag-mag2, 'r*-');
ylabel('幅值差(db)');
legend('ESO', '无ES0');

subplot 212
hold on
grid on
semilogx(data.fre, data.phi-phi1, 'b*-');
semilogx(data.fre, data.phi-phi2, 'r*-');
ylabel('幅值(db)');
xlabel('角频率(rad/s)');








