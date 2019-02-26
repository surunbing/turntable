% clc, clear
% close all

% 添加路径
% projectPath = pwd;
% addpath(genpath(projectPath)); % Add project folder and subfolders to path
% rmpath(genpath([projectPath,'/.git/'])); % remove git from matlab path
% savepath;
% 
% parameter_init;
% global parameter

bforward = 0;

ratio = parameter.ratio;
bandwidth = parameter.bandwidth;
phi_advance = parameter.phi_advance;
phi_advance_margin = parameter.phi_advance_margin;

%% Object
K = parameter.K;%1.56 * 180 / pi;
taue = parameter.taue;%0.00396;% 0.0039035;
taum = parameter.taum;%0.09947;%0.984871194396488;
G = tf(K, [taue * taum, taum + taue, 1, 0]);
%% 加入惯性环节
T = 1 / (bandwidth / 2 / pi) / parameter.Tratio; 
Inertial = tf(1, [T, 1]);

%% 给定频率点
fre = logspace(-1, 2.3) * 2 * pi;  %% 1 - 100
[mag, phi] = bode(G * Inertial, fre);

data.fre = fre;
data.mag = zeros(length(fre), 1);
data.phi = zeros(length(fre), 1);

for i = 1:length(fre)
    data.mag(i) = 20 * log10(mag(1, 1, i));
    data.phi(i) = phi(1, 1, i);
end

%% 补相位

advance = FillPhase(data, ratio, phi_advance, phi_advance_margin, bandwidth);

P = advance.P * Inertial;

figurename('前向通道');
margin(P * G);
grid on;
figurename('前向通道闭环');
bode(P * G / (1 + P * G));
grid on

Design_Lowgain;

phi_creg = parameter.phi_creg;
mag_creg = parameter.mag_creg;
flag_add = 1; % 1: mag, 2: phi
bfailure_pre = 0;
trap_pre = 0;
later_pre = 0;
%% 是否需要加入能否设计出的评估
while 1
    [trap, later, bfailure, data_check, num] = wave_repair(P, G, 0, phi_creg, mag_creg);
    if bfailure == -1 && bfailure_pre == 1
        break;
    end
    bfailure_pre = bfailure;
    if bfailure == 1
        trap_pre = trap;
        later_pre = later;
        %% 性能调优，力图做到更好的数值指标
        if num <= 2
            num_max = min(3, num + 1);
        end
        if flag_add == 1
            mag_creg = mag_creg / parameter.rdiv;
            flag_add = 2;
        elseif flag_add == 2
            phi_creg = phi_creg / parameter.rdiv;
            flag_add = 1;
        end
    else 
        %% 放宽要求，为了前馈和前置滤波器做准备，其中相频优先保证
         if flag_add == 1
            mag_creg = mag_creg * parameter.rdiv;
            flag_add = 2;
        elseif flag_add == 2
            phi_creg = phi_creg * parameter.rdiv;
            flag_add = 1;
        end
    end
end


P = P * later_pre.G * Glow;
for i = 1 : trap_pre.num
    P = P * trap_pre.G(i);
end
figurename('陷波滤波器');
margin(P * G);
grid on
hold on
margin(P * G / Glow);

figurename('陷波滤波器闭环');
bode(P * G / (1 + P * G));
grid on

K_model = parameter.K;
taum = parameter.taum;
taue = parameter.taue;
bforward = 0;
if mag_creg > parameter.maglim  || phi_creg > parameter.philim
    bforward = 1;
    option.type = 'transfer-function';
    [forward, exitflag] = design_forward(P, G, option);
    figurename('顺馈');
    bode((P * G + G * forward.G)/ (1 + P * G));
    grid on 
end

TSp = 0.0005;
fid = fopen(outputfile, 'wt+');
%% 指令预处理
fprintf(fid, '[DOF_%d]\n', nDOF);
fprintf(fid, '// %d-%d 为指令预处理环节，不用时，第1个参数为，其它为0 \n', nZLYCLStart, nZLYCLEnd);
for i = nZLYCLStart : 1 : nZLYCLEnd
    fprintf(fid, 'Link_%02d=%.12f, %.12f, %.12f, %.12f, %.12f, \n', i, 1.0, 0.0, 0.0, 0.0, 0.0);
end

fprintf(fid, '// %d-%d 为顺馈环节，link%d为增益 不用时全部参数必须为0\n', nQKStart, nQKEnd, nQKEnd);
if bforward == 1
    a = [taum, 1];
    b = [taue, 1];
    a_aux = [1 / (parameter.para_aux1 * 2 * pi), 1];
    b_aux = [1 / (parameter.para_aux2 * 2 * pi), 1];
    [dNumf, dDenf] = c2dm(a, a_aux, TSp, 'tustin');
    fprintf(fid, 'Link_%02d=%.l2f,%.12f,%.12f, %.12f,%.12f,\n', nQKStart, -dDenf(2),dNumf(1),dNumf(2), 0, 0);
    [dNumf, dDenf] = c2dm(b, b_aux, TSp, 'tustin');
    fprintf(fid, 'Link_%02d=%.12f, %.12f, %.12f, %.12f, %.12f,\n', nQKStart + 1, -dDenf(2),dNumf(1),dNumf(2), 0, 0);
    fprintf(fid, 'Link_%02d=%.12f, %.12f, %.12f, %.12f, %.12f,\n', nQKEnd, forward.K / K_model, 0, 0, 0, 0);
else
    for i = nQKStart : 1 : nQKEnd
        fprintf(fid, 'Link_%02d=%.12f, %.12f, %.12f, %.12f, %.12f,\n', i, 1.0, 0.0, 0.0, 0.0, 0.0);
    end
end

fprintf(fid, '// %d-%d 为串联校正环节，link%d为增益 不用时全部参数必须为0\n', nJZStart, nJZEnd, nJZEnd);
num = 0;

%% 输出controller
% 惯性环节
[dNumd,dDend] = c2dm(1, Inertial.Denominator{1,1}, TSp, 'tustin');
fprintf(fid, 'Link_%02d=%.12f, %.12f, %.12f, %.12f, %.12f,\n', nJZStart + num, -dDend(2), dNumd(1),dNumd(2), 0, 0);
num = num + 1;
%补充相位
for i = 1 : advance.count + 1
    a = advance.G(i).Numerator{1, 1};
    b = advance.G(i).Denominator{1, 1};
    [dNuml,dDenl] = c2dm(a, b, TSp, 'tustin');
    fprintf(fid, 'Link_%02d=%.12f, %.12f, %.12f, %.12f, %.12f,\n', nJZStart + num + i- 1, -dDenl(2),dNuml(1),dNuml(2), 0, 0);
end
num = num + advance.count + 1;

% later
a = later_pre.G.Numerator{1, 1};
b = later_pre.G.Denominator{1, 1};
[dNuml,dDenl] = c2dm(a, b, TSp, 'tustin');
fprintf(fid, 'Link_%02d=%.12f, %.12f, %.12f, %.12f, %.12f,\n', nJZStart + num, -dDenl(2),dNuml(1),dNuml(2), 0, 0);
num = num + 1;
% trap
for i = 1 : trap_pre.num
    a = [1, trap_pre.e(i) * trap_pre.T(i), trap_pre.f(i) * trap_pre.f(i)];
    b = [1, trap_pre.T(i), trap_pre.f(i) * trap_pre.f(i)];
    [dNumt,dDent] = c2dm(a, b, TSp, 'tustin');
    fprintf(fid, 'Link_%02d=%.12f, %.12f, %.12f, %.12f, %.12f,\n', nJZStart + num + i - 1, -dDent(2),-dDent(3),dNumt(1),dNumt(2),dNumt(3));
end
num = num + trap_pre.num;

%% 迟后
for i = 1 : LowGain.count
    tau = 1 / (sqrt(LowGain.alpha(i)) * LowGain.fre(i));
    [dNuml,dDenl] = c2dm([tau, 1], [LowGain.alpha(i) * tau, 1], TSp, 'tustin');
    fprintf(fid, 'Link_%02d=%.12f, %.12f, %.12f, %.12f, %.12f,\n', nJZStart + num + i - 1, -dDenl(2),dNuml(1),dNuml(2), 0, 0);
end
num = num + LowGain.count;
fprintf(fid, 'Link_%02d=%.12f, %.12f, %.12f, %.12f, %.12f,\n',nJZStart + num, LowGain.K, 0, 0, 0, 0);
num = num + 1;
fprintf(fid, 'Link_%02d=%.12f, %.12f, %.12f, %.12f, %.12f,\n',nJZStart + num, advance.gain, 0, 0, 0, 0);
num = num + 1;
for i = num + nJZStart : 1 : nJZEnd
    fprintf(fid, 'Link_%02d=%.12f, %.12f, %.12f, %.12f, %.12f,\n',i, 1, 0, 0, 0, 0);
end
fclose(fid);

autoArrangeFigures



