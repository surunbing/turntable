clc, clear
close all

fre_start = 100;
fre_end = 150;

fre_array = fre_start : 1 : fre_end;

turntable_bode.fre = fre_array * 2 * pi;
turntable_bode.mag = zeros(length(fre_array), 1);
turntable_bode.phi = zeros(length(fre_array), 1);
Type = 1; %% sine
bTf = 1;  %% ��Ħ��
i = 1;
%% Sweep
for fre = fre_array
    
    t_off = 10;
    mag = 5;
    fre = fre_array(i);
    Tsim = t_off;
    sim('Turntable_sweep.slx',Tsim);

    figure(1)
    subplot 211
    plot(in_out(2000:18000, 1), in_out(2000:18000, 2), 'r-');
    hold on
    grid on
    plot(in_out(2000:18000, 1), in_out(2000:18000, 3), 'b-');
%     
    % FFT
    y = fft(in_out(2000:18000, 2));
    y1 = fft(in_out(2000:18000, 3));
    fs = 2000;       %%����Ƶ��
    n =0: 1: length(in_out) - 1;
    N = length(in_out) - 4000;
    f = n * fs / N;    %Ƶ������
    Mag=abs(y);
    Mag1 = abs(y1);
    subplot 212
    plot(f(1: fix(N / 2)),Mag(1: fix(N / 2)) * 2 / N, 'r-');
    hold on
    grid on
    plot(f(1: fix(N / 2)),Mag1(1: fix(N / 2)) * 2 / N, 'b-');
    
    num = fre * N / (n * fs);
    
    turntable_bode.mag(i) = 20 * (log10(max(Mag1(2 : end))) - log10(max(Mag(2 : end))));
    num = find(Mag1 == max(Mag1(2 : end)));
    num = num(1);
    ang = angle(y1(num));
    if(ang < 0)
        ang = ang + 2 * pi;
    end
    ang1 = angle(y(num));
    if(ang1 < 0)
        ang1 = ang1 + 2 * pi;
    end
    turntable_bode.phi(i) = (ang - ang1) / pi * 180;
    i = i + 1;    
%     close all;
    i;
end
figure(3);
semilogx(turntable_bode.fre, turntable_bode.mag, 'r*-');
grid on
figure(2);
semilogx(turntable_bode.fre, turntable_bode.phi, 'r*-');
grid on

%% ��ƶ���
% figure(4);
% K = 1.56 * 180 / pi;
% taue = 0.0039035;
% taum = 0.984871194396488;
% G = tf(K, [taue * taum, taum, 1, 0]);
% bode(G);
% grid on
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