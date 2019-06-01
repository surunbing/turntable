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

dataG.fre = linspace(1, 70, 70) * 2 * pi;
[mag, phi] = bode(G, dataG.fre);
dataG.mag = 20 * log10(reshape(mag, [length(dataG.fre), 1]));
dataG.phi = reshape(phi, [length(dataG.fre), 1]);

k_search = 0;
mag_creg_up = mag_creg;
mag_creg_lb = 0.001;
phi_creg_up = phi_creg;
phi_creg_lb = 0.001;

mag_creg_pre = 0;
phi_creg_pre = 0;
bSearch_up = 1;

while 1
    [trap, later, G_Inertial, bfailure, data_check, num] = wave_repair(P * Glow, dataG, 0, phi_creg, mag_creg);
    
    
    k_search = [k_search, phi_creg];
    mag_creg_pre = mag_creg;
    phi_creg_pre = phi_creg;
    trap_pre = trap;
    later_pre = later;
    %% 找可以的上限 二分法边界
    if bfailure ~= 1 && bSearch_up == 1
        mag_creg = mag_creg * 1.5;
        phi_creg = phi_creg * 1.5;
    elseif bfailure == 1 && bSearch_up == 1
        bSearch_up = 0; %% 转入二分查找
        mag_creg_up = mag_creg;
        phi_creg_up = phi_creg;
        mag_creg = (mag_creg_up + mag_creg_lb) / 2;
        phi_creg = (phi_creg_up + phi_creg_lb) / 2;
        continue;
    end
    if bSearch_up == 0 && bfailure ~= 1
        mag_creg_lb = mag_creg;
        phi_creg_lb = phi_creg;
        mag_creg = (mag_creg_up + mag_creg_lb) / 2;
        phi_creg = (phi_creg_up + phi_creg_lb) / 2;
    elseif bSearch_up == 0 && bfailure == 1
        mag_creg_up = mag_creg;
        phi_creg_up = phi_creg;
        mag_creg = (mag_creg_up + mag_creg_lb) / 2;
        phi_creg = (phi_creg_up + phi_creg_lb) / 2;
        if abs(mag_creg_pre - mag_creg) < 0.1 && abs(phi_creg_pre - phi_creg) < 0.1
            break;
        end
    end
%     
%     if bfailure == -1 && bfailure_pre == 1
%         break;
%     end
%     bfailure_pre = bfailure;
%     if bfailure == 1
%         trap_pre = trap;
%         later_pre = later;
%         %% 性能调优，力图做到更好的数值指标
%         if num <= 2
%             num_max = min(3, num + 1);
%         end
%         if flag_add == 1
%             mag_creg = mag_creg / parameter.rdiv;
%             flag_add = 2;
%         elseif flag_add == 2
%             phi_creg = phi_creg / parameter.rdiv;
%             flag_add = 1;
%         end
%     else 
%         %% 放宽要求，为了前馈和前置滤波器做准备，其中相频优先保证
%          if flag_add == 1
%             mag_creg = mag_creg * parameter.rdiv;
%             flag_add = 2;
%         elseif flag_add == 2
%             phi_creg = phi_creg * parameter.rdiv;
%             flag_add = 1;
%         end
%     end
end
K2 = P * G_Inertial * later_pre.G;
K = P * later_pre.G * G_Inertial;
for i = 1 : trap_pre.num
    K = K * trap_pre.G(i);
end
figurename('陷波滤波器2');
margin(K2* G);
grid on
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
    [forward, exitflag] = design_forward(K, option);
%     [forward, exitflag] = design_forward(K, G, option);
    figurename('顺馈');
    G1 = tf(K_model, [taum * taue, taum + taue, 1, 0]);
    bode((K * G1 + G1 * forward.G)/ (1 + K * G1));
    grid on
end


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

%% 迟后
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
