function [trap, later, bfailure, data_check] = wave_repair(P, G, para, wc_up, data, bandwidth)
%   repair 波形，通过调整wc， 期望角的裕度， trap个数 与约束范围来调整波形到双十指标内
ratio_max = wc_up / bandwidth;
%% 是否需要加入能否设计出的评估
ratio = wc_up / bandwidth * 0.8;
bandwidth1 = max([bandwidth + pi, para.dt]);
phi_margin = 122;
phi_reg = 8;
num = 1;
num_max = 4;

flag_add = 1; % 1: ratio, 2: phi_reg
bfailure = 1;
while 1
   [later, fval, exitflag] = Holddonewc(P, G, para, bandwidth1, bandwidth * ratio, phi_margin);
   [Gm, Pm, Wgm, Wpm] = margin(P * G * later.G);
   phi_diff = Pm - 45;
   phi_reg = min(phi_diff, phi_reg);      
   [trap, fval, exitflag] = trapdesign(P * later. G, G, bandwidth, num, phi_reg, Wpm);
   if exitflag == 1 || exitflag == 2
       bfailure = 1;
       break;
       %% 检查效果，找到更好的解
        %% 维持现在的参数， 修改约束使得性能更好
   elseif exitflag == 0 || exitflag == -2
       if (ratio < ratio_max || phi_reg < phi_diff) && num <= num_max
           if flag_add == 1
               ratio = min(ratio * 1.025, ratio_max);
               flag_add = 2;
           elseif flag_add == 2
               phi_reg = min(phi_reg + 0.5, phi_diff);
               flag_add = 1;
           end
       elseif (ratio >= ratio_max && phi_reg >= phi_diff) && num <= num_max
           num = min(num + 1, num_max);
           ratio = wc_up / bandwidth * 0.85;
           phi_reg = 8;
           flag_add = 1;
       else 
           bfailure = -1;
           break;
       end      
   end   
end
data_check = 0;
%       K = P * G * later.G;
%    for i = 1 : num
%         K = K * trap.G(i); 
%    end



% figurename('迟后');
% margin(P * G * later.G);
% grid on
% figurename('迟后闭环');
% bode(P * G * later.G / (1 + P * G * later.G));
% grid on

% num = 4;
% K = P * G * later.G;
% for i = 1 : num
%    K = K * trap.G(i); 
% end
% figurename('陷波滤波器');
% margin(K);
% grid on
% figurename('陷波滤波器闭环');
% bode(K / (1 + K));
% grid on
end

