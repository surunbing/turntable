function [trap, later, G_Inertial, bfailure, data_check, num] = wave_repair(P, G, para, phi_creg, mag_creg)
%   repair 波形，通过调整wc， 期望角的裕度， trap个数 与约束范围来调整波形到双十指标内
global parameter bflag_tradition
bandwidth = parameter.bandwidth;
ratio_max = parameter.ratio;
%% 是否需要加入能否设计出的评估
ratio = parameter.start_ratio;
bandwidth1 = parameter.bandwidth + pi; %max([parameter.bandwidth + pi, para.dt]);
phi_margin = parameter.phi_margin;
phi_reg = parameter.phi_reg;
num_max = parameter.num_max;
num = 1;

if bflag_tradition == 1
    G_Inertial = 1;
else
    % 增加一个惯性环节进行低通滤波
    % 增加一个超前环节改善惯性环节引起的相位裕度不足
    fre_Inertial = parameter.Tratio * parameter.bandwidth;
    G_Inertial = tf(1, [1 / fre_Inertial, 1]);
    % 计算损失的相位裕度
    [~, phi_Inertial_cost] = bode(G_Inertial, bandwidth * ratio);
    %计算补充相角的超前环节
    phi_Inertial_cost = abs(phi_Inertial_cost);
    m = sin((phi_Inertial_cost + 3) / 180 * pi);
    alpha = (1 - m) / (1 + m);
    tau = 1 / (sqrt(alpha) * bandwidth * ratio);
    T = alpha * tau;
    G_advance = tf([tau, 1], [T, 1]);
    % 计算增益
    G_Inertial = G_Inertial * G_advance * sqrt(alpha);
end

flag_add = 1; % 1: ratio, 2: phi_reg
bfailure = 0;
later_pre = 1;
while 1
   [later, fval, exitflag] = Holddonewc(P, G, bandwidth1, bandwidth * ratio, phi_margin);
   if exitflag ~= -2
       later_pre = later;
   else
       later = later_pre;
   end
   [Gm, Pm, Wgm, Wpm] = margin(P * G * later.G * G_Inertial);
   phi_diff = Pm - parameter.phimarginmin;
   phi_reg = min(phi_diff, phi_reg);   
   [~, phi_marginreg] = bode(P * G * later.G * G_Inertial, bandwidth);
   phi_marginreg = min(42, 177 - abs(phi_marginreg)); %  需要修改
   [trap, fval, exitflag] = trapdesign(P * later.G * G_Inertial, G, bandwidth, num, phi_reg, Wpm, phi_creg, mag_creg, phi_marginreg);
   if exitflag == 1 || exitflag == 2
       bfailure = 1;
       break;
       %% 检查效果，找到更好的解  性能调优： 放在外层的循环
        %% 维持现在的参数， 修改约束使得性能更好
   elseif exitflag == 0 || exitflag == -2
       if (ratio < ratio_max || phi_reg < phi_diff) && num <= num_max
           if flag_add == 1
               ratio = min(ratio * parameter.ratiodiv, ratio_max);
               flag_add = 2;
           elseif flag_add == 2
               phi_reg = min(phi_reg + parameter.phidiv, phi_diff);
               flag_add = 1;
           end
       elseif (ratio >= ratio_max && phi_reg >= phi_diff) && num <= num_max
           num = num + 1; %min(num + 1, num_max);
           ratio = parameter.start_ratio;
           phi_reg = parameter.phi_reg;
           flag_add = 1;
       end      
       if num > num_max
           num = num_max;
           bfailure = -1;

           break;
       end
   end   
end
data_check = 0;
end

