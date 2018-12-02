% close all
% 
% parameter_init;
% 
% global parameter


[P, G, para] = direct_design();

figurename('直接设计开环');
margin(P * G);
grid on

figurename('直接设计闭环');
bode(P * G / (1 + P * G));
grid on

Design_Lowgain

%% 衡量约束的量
phi_creg = parameter.phi_creg;
mag_creg = parameter.mag_creg;
flag_add = 1; % 1: mag, 2: phi
bfailure_pre = 0;
trap_pre = 0;
later_pre = 0;
%% 是否需要加入能否设计出的评估
while 1
    [trap, later, bfailure, data_check, num] = wave_repair(P, G, para, phi_creg, mag_creg);
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

K = P * later_pre.G * Glow;
for i = 1 : trap_pre.num
    K = K * trap_pre.G(i);
end
figurename('陷波滤波器');
margin(K* G);
grid on
figurename('陷波滤波器闭环');
bode(K * G / (1 + K * G));
grid on

%% 顺馈
K_model = parameter.K;
taum = parameter.taum;
taue = parameter.taue;
bforward = 0;
if mag_creg > parameter.maglim  || phi_creg > parameter.philim
    bforward = 1;
    option.type = 'transfer-function';
    [forward, exitflag] = design_forward(K, G, option);
    figurename('顺馈');
    G1 = tf(K_model, [taum * taue, taum + taue, 1, 0]);
    bode((K * G + G * forward.G)/ (1 + K * G));
    grid on
end


fid = fopen('C:\Users\Momenta\Documents\毕业设计\CSDA_FANGXUN\turntable\controller.txt', 'wt+');
%% 指令预处理
fprintf(fid, '[DOF_%d]', nDOF);
fprintf(fid, '// %d-%d 为指令预处理环节，不用时，第1个参数为，其它为0', nZLYCLStart, nZLYCLEnd);
for i = nZLYCLStart : 1 : nZLYCLEnd
    fprintf(fid, 'Link_%02d=%lf, %lf, %lf, %lf, %lf,\n', i, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
end

fprintf(fid, '// %d-%d 为顺馈环节，link6为增益 不用时全部参数必须为0', nQKStart, nQKEnd);
if bforward == 1
    a = [taum, 1];
    b = [taue, 1];
    a_aux = [1 / (parameter.para_aux1 * 2 * pi), 1];
    b_aux = [1 / (parameter.para_aux2 * 2 * pi), 1];
    [dNumf, dDenf] = c2dm(a, a_aux, TSp, 'tustin');
    fprintf(fid, 'Link_%02d=%lf, %lf, %lf, %lf, %lf,\n', nQKStart, -dDenf(2),dNumf(1),dNumf(2), 0, 0);
    fprintf(fid, '%.12f, %.12f, %.12f, %.12f, %.12f, 前馈\n', -dDenf(2),dNumf(1),dNumf(2), 0, 0);
    [dNumf, dDenf] = c2dm(b, b_aux, TSp, 'tustin');
    fprintf(fid, '%.12f, %.12f, %.12f, %.12f, %.12f, 前馈\n', -dDenf(2),dNumf(1),dNumf(2), 0, 0);
    fprintf(fid, '%.12f, %.12f, %.12f, %.12f, %.12f, 前馈增益\n',forward.K / K_model,0,0,0,0);
else
    for i = nQKStart : 1 : nQKEnd
        fprintf(fid, 'Link_%02d=%lf, %lf, %lf, %lf, %lf,\n', i, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    end
end


T = para.T;
omegan = para.omegan;
xi = para.xi;
a = conv([taue, 1], [taum, 1]);
b = [T, 2 * T * xi * omegan + 1, omegan * (T * omegan + 2 * xi)];
k = omegan * omegan / K_model;
TSp = 0.0005;
[dNumd,dDend] = c2dm(a, b, TSp, 'tustin');
fprintf(fid, '%.12f, %.12f, %.12f, %.12f, %.12f, 直接增益\n',k,0,0,0,0);
fprintf(fid, '%.12f, %.12f, %.12f, %.12f, %.12f, 直接\n',-dDend(2),-dDend(3),dNumd(1),dNumd(2),dNumd(3));

% later
a = later_pre.G.Numerator{1, 1};
b = later_pre.G.Denominator{1, 1};
[dNuml,dDenl] = c2dm(a, b, TSp, 'tustin');
fprintf(fid, '%.12f, %.12f, %.12f, %.12f, %.12f, 迟后\n', -dDenl(2),dNuml(1),dNuml(2), 0, 0);
% trap
for i = 1 : trap_pre.num
    a = [1, trap_pre.e(i) * trap_pre.T(i), trap_pre.f(i) * trap_pre.f(i)];
    b = [1, trap_pre.T(i), trap_pre.f(i) * trap_pre.f(i)];
    [dNumt,dDent] = c2dm(a, b, TSp, 'tustin');
    fprintf(fid, '%.12f, %.12f, %.12f, %.12f, %.12f, 陷波\n',-dDent(2),-dDent(3),dNumt(1),dNumt(2),dNumt(3));
end

%% 迟后
for i = 1 : LowGain.count
    tau = 1 / (sqrt(LowGain.alpha(i)) * LowGain.fre(i));
    [dNuml,dDenl] = c2dm([tau, 1], [LowGain.alpha(i) * tau, 1], TSp, 'tustin');
    fprintf(fid, '%.12f, %.12f, %.12f, %.12f, %.12f, 低频迟后\n',-dDenl(2),dNuml(1),dNuml(2), 0, 0);
end
fprintf(fid, '%.12f, %.12f, %.12f, %.12f, %.12f, 低频增益\n',LowGain.K,0,0,0,0);

fclose(fid);
autoArrangeFigures;
