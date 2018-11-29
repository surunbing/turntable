clc, clear
close all
nlength = 75;
frequence = linspace(1, nlength, nlength) * 2 * pi;
ET205 = load('ET205主轴3Vsweepa.csv');
mag = ET205(1: nlength, 1);
phi = ET205(1: nlength, 2);
response = mag.*exp(1j*phi*pi/180);
fr_data = idfrd(response,frequence,0);
opt = procestOptions('Focus','simulation','Display','on','SearchMethod', 'fmincon');
sysP1D_noise = procest(fr_data,'P2I',opt);
G = tf(sysP1D_noise.Kp, [sysP1D_noise.Tp1 * sysP1D_noise.Tp2, sysP1D_noise.Tp1 + sysP1D_noise.Tp2, 1, 0]);
Mag_real = 20 * log10(mag);
Phi_real = phi;
Mag_fit = zeros(length(frequence), 1);
Phi_fit = zeros(length(frequence), 1);
[mag, phi] = bode(G, frequence);
for i = 1 : length(frequence)
    Mag_fit(i) = 20 * log10(mag(1, 1, i));
    Phi_fit(i) = phi(i);
end
figurename('对象对比2');
subplot 211;
semilogx(frequence, Mag_real, 'r*-');
hold on
semilogx(frequence, Mag_fit, 'g*-');
grid on

subplot 212;
semilogx(frequence, Phi_real, 'r*-');
hold on
semilogx(frequence, Phi_fit, 'g*-');
grid on

autoArrangeFigures;