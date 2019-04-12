bforward = 0;

ratio = parameter.ratio;
bandwidth = parameter.bandwidth;
phi_advance = parameter.phi_advance;
phi_advance_margin = parameter.phi_advance_margin;

%% Object
%% 加载数据
% ET205 = load('ET205主轴3Vsweep-test.csv');
G_model = tf(K, [taue * taum taue + taum 1 0]);
data.fre = linspace(1, 70, 70)' * 2 * pi;
[mag, phi] = bode(G_model, data.fre);
data.mag = 20 .* log10(reshape(mag, [70, 1]));
data.phi = reshape(phi, [70, 1]);

dataG = data;
complex_G = 10 .^ (data.mag / 20) .* complex(cos(data.phi / 180  * pi), sin(data.phi / 180  * pi));


figurename('对象特性');
subplot 211 
semilogx(dataG.fre, dataG.mag, 'r*-');
grid on
subplot 212
semilogx(dataG.fre, dataG.phi, 'r*-');
grid on


G = tf(K, [taue * taum, taum + taue, 1, 0]);
%% 加入惯性环节
T = 1 / (bandwidth / 2 / pi) / parameter.Tratio; 
Inertial = tf(1, [T, 1]);

%% 给定频率点
[mag, phi] = bode(Inertial, data.fre);
mag = 20 * log10(reshape(mag, [length(data.fre), 1]));
phi = reshape(phi, [length(data.fre), 1]);
data.mag = mag + data.mag;
data.phi = phi + data.phi;

%% 补相位
advance = FillPhase(data, ratio, phi_advance, phi_advance_margin, bandwidth);

P = advance.P * Inertial;

figurename('前向通道');
[mag, phi] = bode(P, data.fre);
mag = 20 * log10(reshape(mag, [length(data.fre), 1]));
phi = reshape(phi, [length(data.fre), 1]);
subplot 211
semilogx(dataG.fre, dataG.mag + mag, 'r-');
grid on;
subplot 212
semilogx(dataG.fre, dataG.phi + phi, 'r-');
grid on
Design_Lowgain;

phi_creg = parameter.phi_creg;
mag_creg = parameter.mag_creg;
trap_pre = 0;
later_pre = 0;
bSearch_up = 1;

mag_creg_up = mag_creg;
mag_creg_lb = 0.001;
phi_creg_up = phi_creg;
phi_creg_lb = 0.001;

mag_creg_pre = 0;
phi_creg_pre = 0;

k_search = 0;
%% 是否需要加入能否设计出的评估
while 1
    [trap, later, ~, bfailure, data_check, num] = wave_repair(P * Glow, dataG, 0, phi_creg, mag_creg);
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
end


P = P * later_pre.G * Glow;
for i = 1 : trap_pre.num
    P = P * trap_pre.G(i);
end
figurename('陷波滤波器');
data_out = translate_data(dataG, P);
subplot 211
semilogx(data_out.fre, data_out.mag, 'r-');
grid on
% margin(P * G);
subplot 212
semilogx(data_out.fre, data_out.phi, 'r-');
grid on

figurename('陷波滤波器闭环');
%% 计算闭环
complex_open = 10 .^ (data_out.mag / 20) .* complex(cos(data_out.phi / 180  * pi), sin(data_out.phi / 180  * pi));
subplot 211
mag = 20 * log10(abs(complex_open ./ (1 + complex_open)));
semilogx(data_out.fre, mag, 'r-');
grid on
subplot 212
phi = angle(complex_open ./ (1 + complex_open)) / pi * 180;
for i = 1 : length(phi)
    if phi(i) > 30
        phi(i) = phi(i) - 360;
    end
end
semilogx(data_out.fre, phi, 'r-');
grid on

K_model = parameter.K;
taum = parameter.taum;
taue = parameter.taue;
bforward = 0;
G = tf(K, [taue * taum, taue + taum, 1, 0]);
if mag_creg > parameter.maglim  || phi_creg > parameter.philim
    bforward = 1;
    option.type = 'transfer-function';
    [forward, exitflag] = design_forward(P, G, option);
    figurename('顺馈');
    bode((P * G + G * forward.G)/ (1 + P * G));
    grid on 
    %% 绘出离散的图
    [mag, phi] = bode(forward.G, data_out.fre);
    mag = reshape(mag, [length(data_out.fre), 1]);
    phi = reshape(phi, [length(data_out.fre), 1]);
    complex_forward = mag .* complex(cos(phi / 180 * pi), sin(phi / 180 * pi));
    complex_close = (complex_open + complex_G .* complex_forward) ./ (1 + complex_open);
    figurename('顺馈离散');
    subplot 211
    semilogx(data_out.fre, 20 * log10(abs(complex_close)), 'r-');
    grid on
    subplot 212
    semilogx(data_out.fre, angle(complex_close) / pi * 180, 'r-');
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



