clc, clear
close all

fre_start = 1;
fre_end = 50;

fre_array = fre_start : 2 : fre_end;

turntable_bode.fre = fre_array * 2 * pi;
turntable_bode.mag = zeros(length(fre_array), 1);
turntable_bode.phi = zeros(length(fre_array), 1);
turntable_bode.magc = zeros(length(fre_array), 1);
turntable_bode.phic = zeros(length(fre_array), 1);
Type = 1; %% sine
bTf = 1;  %% 0无摩擦  1有摩擦
bFt = 0;  %% 0 有力矩波动  1 无力矩波动
i = 1;
%% Sweep
for fre = fre_array
    
    t_off = 10;
    mag = 5;
    fre = fre_array(i);
    Tsim = t_off;
    sim('Turntable_sweep.slx',Tsim);

    turntable_bode.magc(i) = 20 * log10(mag);
    turntable_bode.phic(i) = 0;
%     figure(1)
%     subplot 211
%     plot(in_out(2000:18000, 1), in_out(2000:18000, 2), 'r-');
%     hold on
%     grid on
%     plot(in_out(2000:18000, 1), in_out(2000:18000, 3), 'b-');
%     
    % FFT
    y = fft(in_out(2000:18000, 2));
    y1 = fft(in_out(2000:18000, 3));
    fs = 2000;       %%采样频率
    n =0: 1: length(in_out) - 1;
    N = length(in_out) - 4000;
    f = n * fs / N;    %频率序列
    Mag=abs(y);
    Mag1 = abs(y1);
%     subplot 212
%     plot(f(1: fix(N / 2)),Mag(1: fix(N / 2)) * 2 / N, 'r-');
%     hold on
%     grid on
%     plot(f(1: fix(N / 2)),Mag1(1: fix(N / 2)) * 2 / N, 'b-');
    
%     num = fre * N / (n * fs);
    n_f = find(Mag == max(Mag));
    n_f = n_f(1);
%     n_f = fre * N / fs;
    turntable_bode.mag(i) = 20 * (log10(Mag1(n_f)) - log10(Mag(n_f)));
%     turntable_bode.mag(i) = 20 * (log10(max(Mag1(2 : end))) - log10(max(Mag(2 : end))));
%     num = find(Mag1 == max(Mag1(2 : end)));
%     num = num(1);
    ang = angle(y1(n_f));
    if(ang < 0)
        ang = ang + 2 * pi;
    end
    ang1 = angle(y(n_f));
    if(ang1 < 0)
        ang1 = ang1 + 2 * pi;
    end
    turntable_bode.phi(i) = (ang - ang1) / pi * 180;
    i = i + 1;    
%     close all;
    i;
end
h1 = figurename('mag');
semilogx(turntable_bode.fre, turntable_bode.mag, 'r*-');
hold on
grid on
h2 = figurename('phi');
semilogx(turntable_bode.fre, turntable_bode.phi, 'r*-');
hold on
grid on


%% 标称对象
K = 1.56 * 180 / pi;
taue = 0.0039035;
taum = 0.984871194396488;
G = tf(K, [taue * taum, taum, 1, 0]);
[mag, phi] = bode(G, turntable_bode.fre);
Mag = zeros(length(turntable_bode.fre), 1);
Phi = zeros(length(turntable_bode.fre), 1);
for i = 1 : length(turntable_bode.fre)
    Mag(i) = 20 * log10(mag(1, 1, i));
    Phi(i) = phi(1, 1, i);
end
figure(h1);
semilogx(turntable_bode.fre, Mag, 'bo-');
grid on
hold on
figure(h2);
semilogx(turntable_bode.fre, Phi, 'bo-');
grid on
hold on
mag = 10 .^ (turntable_bode.mag / 20);
phi = turntable_bode.phi;
response = mag.*exp(1j*phi*pi/180);
fr_data = idfrd(response,turntable_bode.fre,0);
sysP1D_noise = procest(fr_data,'P2I');

G = tf(sysP1D_noise.Kp, [sysP1D_noise.Tp1 * sysP1D_noise.Tp2, sysP1D_noise.Tp1 + sysP1D_noise.Tp2, 1, 0]);
[mag, phi] = bode(G, turntable_bode.fre);
Mag = zeros(length(turntable_bode.fre), 1);
Phi = zeros(length(turntable_bode.fre), 1);
for i = 1 : length(turntable_bode.fre)
    Mag(i) = 20 * log10(mag(1, 1, i));
    Phi(i) = phi(1, 1, i);
end
% Mag = 20 * log10(abs(response));
% Phi = angle(response) / pi * 180;
% for i = 1 : length(Phi)
%    if Phi(i) > 0
%        Phi(i) = Phi(i) - 360;
%    end
%        
% end
figure(h1);
semilogx(turntable_bode.fre, Mag, 'yo-');
grid on
figure(h2);
semilogx(turntable_bode.fre, Phi, 'yo-');
grid on


% 
% figure(5)
% T = G / (1 + G);
% bode(T);
% grid on
% 
% figure(6)
% S = 1 / (1 + G);
% bode(S);
% grid on
% 
% save('turntable_fre_nodisturb.mat', 'turntable_bode');