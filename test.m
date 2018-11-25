clc, clear;
close all

% K = 434;
% taum = 0.67;
% taue = 0.0035;

K = 1.56 * 180 / pi;
taue = 0.0039035;
taum = 0.984871194396488;

bandwidth = 10 * 2 * pi;
wc_max = 650;
[P, G, para] = direct_design(bandwidth, wc_max, K, taum, taue);

figurename('ֱ����ƿ���');
margin(P * G);
grid on

figurename('ֱ����Ʊջ�');
bode(P * G / (1 + P * G));
grid on

Design_Lowgain

wc_max = bandwidth * 2.1;
%% ����Լ������
phi_creg = 9;
mag_creg = 0.8;
num_max = 3;
flag_add = 1; % 1: mag, 2: phi
bfailure_pre = 0;
trap_pre = 0;
later_pre = 0;
%% �Ƿ���Ҫ�����ܷ���Ƴ�������
while 1
    [trap, later, bfailure, data_check, num] = wave_repair(P * Glow, G, para, wc_max, 0, bandwidth, phi_creg, mag_creg, num_max);
    if bfailure == -1 && bfailure_pre == 1
        break;
    end
    bfailure_pre = bfailure;
    if bfailure == 1
        trap_pre = trap;
        later_pre = later;
        %% ���ܵ��ţ���ͼ�������õ���ֵָ��
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
        %% �ſ�Ҫ��Ϊ��ǰ����ǰ���˲�����׼����������Ƶ���ȱ�֤
         if flag_add == 1
            mag_creg = mag_creg * 1.05;
            flag_add = 2;
        elseif flag_add == 2
            phi_creg = phi_creg * 1.05;
            flag_add = 1;
        end
    end
end

K = P * later_pre.G * Glow;
for i = 1 : trap_pre.num
    K = K * trap_pre.G(i);
end
figurename('�ݲ��˲���');
margin(K* G);
grid on
figurename('�ݲ��˲����ջ�');
bode(K * G / (1 + K * G));
grid on

%% ˳��
if mag_creg > 0.8 || phi_creg > 9
    option.type = 'transfer-function';
    [forward, exitflag] = design_forward(K, G, 0, 0, 0, 0, 0.7, bandwidth, option);
    figurename('˳��');
    bode((K * G + G * forward.G)/ (1 + K * G));
    grid on
end

%% ����Ծ����
t = 0 : 0.0005 : 10;
u = ones(length(t), 1) * 3;
out = lsim((K * G + G * forward.G)/ (1 + K * G), u, t);
out1 = lsim((K * G)/ (1 + K * G), u, t);
figurename('��Ծ');
plot(t, u, 'r');
hold on
grid on
plot(t, out, 'b');
hold on
plot(t, out1, 'g');

% a = omegan * omegan * conv([taue, 1], [taum, 1]);
% b = K * [T, 2 * T * xi * omegan + 1, omegan * (T * omegan + 2 * xi)];
T = para.T;
omegan = para.omegan;
xi = para.xi;
a = conv([taue, 1], [taum, 1]);
b = [T, 2 * T * xi * omegan + 1, omegan * (T * omegan + 2 * xi)];
k = omegan * omegan / K_model;
TSp = 0.0005;
[dNumd,dDend] = c2dm(a, b, TSp, 'tustin');

% later
a = later_pre.G.Numerator{1, 1};
b = later_pre.G.Denominator{1, 1};
[dNuml,dDenl] = c2dm(a, b, TSp, 'tustin');

% trap
for i = trap_pre.count
    
end

autoArrangeFigures;
