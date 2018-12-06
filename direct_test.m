% close all
% 
% parameter_init;
% 
% global parameter


[P, G, para] = direct_design();

figurename('ֱ����ƿ���');
margin(P * G);
grid on

figurename('ֱ����Ʊջ�');
bode(P * G / (1 + P * G));
grid on

Design_Lowgain

%% ����Լ������
phi_creg = parameter.phi_creg;
mag_creg = parameter.mag_creg;
flag_add = 1; % 1: mag, 2: phi
bfailure_pre = 0;
trap_pre = 0;
later_pre = 0;
%% �Ƿ���Ҫ�����ܷ���Ƴ�������
while 1
    [trap, later, bfailure, data_check, num] = wave_repair(P, G, para, phi_creg, mag_creg);
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
K_model = parameter.K;
taum = parameter.taum;
taue = parameter.taue;
bforward = 0;
if mag_creg > parameter.maglim  || phi_creg > parameter.philim
    bforward = 1;
    option.type = 'transfer-function';
    [forward, exitflag] = design_forward(K, G, option);
    figurename('˳��');
    G1 = tf(K_model, [taum * taue, taum + taue, 1, 0]);
    bode((K * G + G * forward.G)/ (1 + K * G));
    grid on
end


fid = fopen(outputfile, 'wt+');
%% ָ��Ԥ����
fprintf(fid, '[DOF_%d]\n', nDOF);
fprintf(fid, '// %d-%d Ϊָ��Ԥ�����ڣ�����ʱ����1������Ϊ������Ϊ0 \n', nZLYCLStart, nZLYCLEnd);
for i = nZLYCLStart : 1 : nZLYCLEnd
    fprintf(fid, 'Link_%02d=%.12f, %.12f, %.12f, %.12f, %.12f, \n', i, 1.0, 0.0, 0.0, 0.0, 0.0);
end

fprintf(fid, '// %d-%d Ϊ˳�����ڣ�link%dΪ���� ����ʱȫ����������Ϊ0\n', nQKStart, nQKEnd, nQKEnd);
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

fprintf(fid, '// %d-%d Ϊ����У�����ڣ�link%dΪ���� ����ʱȫ����������Ϊ0\n', nJZStart, nJZEnd, nJZEnd);
num = 0;

T = para.T;
omegan = para.omegan;
xi = para.xi;
a = conv([taue, 1], [taum, 1]);
b = [T, 2 * T * xi * omegan + 1, omegan * (T * omegan + 2 * xi)];
k = omegan * omegan / K_model;
TSp = 0.0005;
[dNumd,dDend] = c2dm(a, b, TSp, 'tustin');
fprintf(fid, 'Link_%02d=%.12f, %.12f, %.12f, %.12f, %.12f,\n', nJZStart + num, -dDend(2),-dDend(3),dNumd(1),dNumd(2),dNumd(3));
num = num + 1;

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

%% �ٺ�
for i = 1 : LowGain.count
    tau = 1 / (sqrt(LowGain.alpha(i)) * LowGain.fre(i));
    [dNuml,dDenl] = c2dm([tau, 1], [LowGain.alpha(i) * tau, 1], TSp, 'tustin');
    fprintf(fid, 'Link_%02d=%.12f, %.12f, %.12f, %.12f, %.12f,\n', nJZStart + num + i - 1, -dDenl(2),dNuml(1),dNuml(2), 0, 0);
end
num = num + LowGain.count;
fprintf(fid, 'Link_%02d=%.12f, %.12f, %.12f, %.12f, %.12f,\n',nJZStart + num, LowGain.K, 0, 0, 0, 0);
num = num + 1;
fprintf(fid, 'Link_%02d=%.12f, %.12f, %.12f, %.12f, %.12f,\n',nJZStart + num, k, 0, 0, 0, 0);
num = num + 1;
for i = num + nJZStart : 1 : nJZEnd
    fprintf(fid, 'Link_%02d=%.12f, %.12f, %.12f, %.12f, %.12f,\n',i, 1, 0, 0, 0, 0);
end


fclose(fid);
autoArrangeFigures;
