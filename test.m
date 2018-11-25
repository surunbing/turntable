clc, clear;
close all

% K = 434;
% taum = 0.67;
% taue = 0.0035;

K = 1.56 * 180 / pi;
taue = 0.0039035;
taum = 0.984871194396488;

bandwidth = 18 * 2 * pi;
wc_max = 650;
[P, G, para] = direct_design(bandwidth, wc_max, K, taum, taue);

figurename('直接设计开环');
margin(P * G);
grid on

figurename('直接设计闭环');
bode(P * G / (1 + P * G));
grid on

Design_Lowgain

wc_max = bandwidth * 2.1;
%% 衡量约束的量
phi_creg = 9;
mag_creg = 0.8;
num_max = 3;
flag_add = 2; % 1: mag, 2: phi
bfailure_pre = 0;
trap_pre = 0;
later_pre = 0;
%% 是否需要加入能否设计出的评估
while 1
    [trap, later, bfailure, data_check, num] = wave_repair(P * Glow, G, para, wc_max, 0, bandwidth, phi_creg, mag_creg, num_max);
    if bfailure == -1 && bfailure_pre == 1
        break;
    end
    bfailure_pre = bfailure;
    if bfailure == 1
        trap_pre = trap;
        later_pre = later;
        %% 性能调优，力图做到更好的数值指标
        if num <= 2
            num_max = min(2, num + 1);
        end
        if flag_add == 1
            mag_creg = mag_creg / 1.05;
            flag_add = 2;
        elseif flag_add == 2
            phi_creg = phi_creg / 1.05;
            flag_add = 1;
        end
    else 
        %% 放宽要求，为了前馈和前置滤波器做准备，其中相频优先保证
         if flag_add == 1
            mag_creg = mag_creg * 1.05;
            flag_add = 2;
        elseif flag_add == 2
            phi_creg = phi_creg * 1.05;
            flag_add = 1;
        end
    end
end

K = P * G * later_pre.G * Glow;
for i = 1 : trap_pre.num
    K = K * trap_pre.G(i);
end
figurename('陷波滤波器');
margin(K);
grid on
figurename('陷波滤波器闭环');
bode(K / (1 + K));
grid on


% Design_Lowgain
% 
% figurename('低频增益');
% K = K * Glow;
% margin(K);
% grid on
% figurename('低频闭环');
% bode(K / (1 + K));
% grid on


autoArrangeFigures;
