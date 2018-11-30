clc, clear
close all

% ���·��
projectPath = pwd;
addpath(genpath(projectPath)); % Add project folder and subfolders to path
rmpath(genpath([projectPath,'/.git/'])); % remove git from matlab path
savepath;

parameter_init;
global parameter

bforward = 0;

ratio = 2.5;

bandwidth = 18 * 2 * pi;

%% Object
K = 493;%1.56 * 180 / pi;
taue = 0.0004;%0.00396;% 0.0039035;
taum = 2.05;%0.09947;%0.984871194396488;
G = tf(K, [taue * taum, taum, 1, 0]);
%% ������Ի���
T = 1 / (bandwidth / 2 / pi) / 15; 
Inertial = tf(1, [T, 1]);
G = G * Inertial;

%% ����Ƶ�ʵ�
fre = logspace(-1, 2.3) * 2 * pi;  %% 1 - 100
[mag, phi] = bode(G, fre);

data.fre = fre;
data.mag = zeros(length(fre), 1);
data.phi = zeros(length(fre), 1);

for i = 1:length(fre)
    data.mag(i) = 20 * log10(mag(1, 1, i));
    data.phi(i) = phi(1, 1, i);
end

%% ����λ
phi_advance = 30;
phi_advance_margin = 0;
advance = FillPhase(data, ratio, phi_advance, phi_advance_margin, bandwidth);

G_P = G * advance.P;

figurename('ǰ��ͨ��');
margin(G_P);
grid on;
figurename('ǰ��ͨ���ջ�');
bode(G_P / (1 + G_P));
grid on

Design_Lowgain;
P = advance.P * Inertial;
K = 493;%1.56 * 180 / pi;
taue = 0.0004;%0.00396;% 0.0039035;
taum = 2.05;%0.09947;%0.984871194396488;
G = tf(K, [taue * taum, taum + taue, 1, 0]);

phi_creg = parameter.phi_creg;
mag_creg = parameter.mag_creg;
flag_add = 1; % 1: mag, 2: phi
bfailure_pre = 0;
trap_pre = 0;
later_pre = 0;
%% �Ƿ���Ҫ�����ܷ���Ƴ�������
while 1
    [trap, later, bfailure, data_check, num] = wave_repair(P, G, 0, phi_creg, mag_creg);
    if bfailure == -1 && bfailure_pre == 1
        break;
    end
    bfailure_pre = bfailure;
    if bfailure == 1
        trap_pre = trap;
        later_pre = later;
        %% ���ܵ��ţ���ͼ�������õ���ֵָ��
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
        %% �ſ�Ҫ��Ϊ��ǰ����ǰ���˲�����׼����������Ƶ���ȱ�֤
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
figurename('�ݲ��˲���');
margin(P * G);
grid on
figurename('�ݲ��˲����ջ�');
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
    figurename('˳��');
    bode((P * G + G * forward.G)/ (1 + P * G));
    grid on 
end


%% ���controller
TSp = 0.0005;
% ���Ի���
[dNumd,dDend] = c2dm(1, Inertial.Denominator{1,1}, TSp, 'tustin');
fid = fopen('controller.txt', 'wt+');
fprintf(fid, '%.12f, %.12f, %.12f, %.12f, %.12f, ����\n',-dDend(2), dNumd(1),dNumd(2), 0, 0);
%������λ


% later
a = later_pre.G.Numerator{1, 1};
b = later_pre.G.Denominator{1, 1};
[dNuml,dDenl] = c2dm(a, b, TSp, 'tustin');
fprintf(fid, '%.12f, %.12f, %.12f, %.12f, %.12f, �ٺ�\n', -dDenl(2),dNuml(1),dNuml(2), 0, 0);
% trap
for i = 1 : trap_pre.num
    a = [1, trap_pre.e(i) * trap_pre.T(i), trap_pre.f(i) * trap_pre.f(i)];
    b = [1, trap_pre.T(i), trap_pre.f(i) * trap_pre.f(i)];
    [dNumt,dDent] = c2dm(a, b, TSp, 'tustin');
    fprintf(fid, '%.12f, %.12f, %.12f, %.12f, %.12f, �ݲ�\n',-dDent(2),-dDent(3),dNumt(1),dNumt(2),dNumt(3));
end
% forward
if bforward == 1
    a = [taum, 1];
    b = [taue, 1];
    a_aux = [1 / (parameter.para_aux1 * 2 * pi), 1];
    b_aux = [1 / (parameter.para_aux2 * 2 * pi), 1];
    [dNumf, dDenf] = c2dm(a, a_aux, TSp, 'tustin');
    fprintf(fid, '%.12f, %.12f, %.12f, %.12f, %.12f, ǰ��\n', -dDenf(2),dNumf(1),dNumf(2), 0, 0);
    [dNumf, dDenf] = c2dm(b, b_aux, TSp, 'tustin');
    fprintf(fid, '%.12f, %.12f, %.12f, %.12f, %.12f, ǰ��\n', -dDenf(2),dNumf(1),dNumf(2), 0, 0);
    fprintf(fid, '%.12f, %.12f, %.12f, %.12f, %.12f, ǰ������\n',forward.K / parameter.K, 0, 0, 0, 0);
end

%% �ٺ�
for i = 1 : LowGain.count
    tau = 1 / (sqrt(LowGain.alpha(i)) * LowGain.fre(i));
    [dNuml,dDenl] = c2dm([tau, 1], [LowGain.alpha(i) * tau, 1], TSp, 'tustin');
    fprintf(fid, '%.12f, %.12f, %.12f, %.12f, %.12f, ��Ƶ�ٺ�\n',-dDenl(2),dNuml(1),dNuml(2), 0, 0);
end
fprintf(fid, '%.12f, %.12f, %.12f, %.12f, %.12f, ��Ƶ����\n',LowGain.K, 0, 0, 0, 0);

fclose(fid);

autoArrangeFigures



