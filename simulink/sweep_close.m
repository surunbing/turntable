clc, clear
close all

% fre_start = 0;
% fre_end = 100;
% fre_step = 0.5;
circle_num = 30;
% 
% fre_array = fre_start: fre_step: fre_end;
% fre_array(1) = 0.1;

fre_array = logspace(-1, 2);

turntable_bode.fre = fre_array * 2 * pi;
turntable_bode.mag = zeros(length(fre_array), 1);
turntable_bode.phi = zeros(length(fre_array), 1);

i = 1;
%% Sweep
for fre = fre_array

    if(fre > 5)
        circle_num = 100;
    end
    
    t_off = 1 / fre * circle_num;
    mag = 1;
%     fre = 5;
    Tsim = t_off;
    sim('Turntable_sweep_close.slx',Tsim);

%     figure(1)
%     subplot 211
%     plot(in_out(:, 1), in_out(:, 2), 'r-');
%     hold on
%     grid on
%     plot(in_out(:, 1), in_out(:, 3), 'b-');
    
    % FFT
    y = fft(in_out(:, 2));
    y1 = fft(in_out(:, 3));
    fs = 2000;       %%采样频率
    n =0: 1: length(in_out) - 1;
    N = length(in_out);
    f = n * fs / N;    %频率序列
    Mag=abs(y);
    Mag1 = abs(y1);
%     subplot 212
%     plot(f(1: N / 2),Mag(1: N / 2) * 2 / N, 'r-');
%     hold on
%     grid on
%     plot(f(1: N / 2),Mag1(1: N / 2) * 2 / N, 'b-');
    
    turntable_bode.mag(i) = 20 * (log10(max(Mag1(2 : end))) - log10(max(Mag(2 : end))));
    num = find(Mag1 == max(Mag1(2 : end)));
    num = num(1);
    turntable_bode.phi(i) = angle(y1(num)) - angle(y(num));
    i = i + 1;    
%     close all;
    i;
end
figure(3);
semilogx(turntable_bode.fre, turntable_bode.mag, 'r*-');
grid on
figure(2);
semilogx(turntable_bode.fre, turntable_bode.phi, 'r*-');


%% 标称对象
figure(4);
K = 1.56 * 180 / pi;
taue = 0.0039035;
taum = 0.984871194396488;
G = tf(K, [taue * taum, taum, 1, 0]);
bode(G);
grid on

figure(5)
T = G / (1 + G);
bode(T);
grid on

figure(6)
S = 1 / (1 + G);
bode(S);
grid on