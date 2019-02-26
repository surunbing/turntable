function [trap, later, bfailure, data_check, num] = wave_repair(P, G, para, phi_creg, mag_creg)
%   repair ���Σ�ͨ������wc�� �����ǵ�ԣ�ȣ� trap���� ��Լ����Χ���������ε�˫ʮָ����
global parameter
bandwidth = parameter.bandwidth;
ratio_max = parameter.ratio;
%% �Ƿ���Ҫ�����ܷ���Ƴ�������
ratio = parameter.start_ratio;
bandwidth1 = parameter.bandwidth + pi; %max([parameter.bandwidth + pi, para.dt]);
phi_margin = parameter.phi_margin;
phi_reg = parameter.phi_reg;
num_max = parameter.num_max;
num = 1;

flag_add = 1; % 1: ratio, 2: phi_reg
bfailure = 0;
while 1
   [later, fval, exitflag] = Holddonewc(P, G, bandwidth1, bandwidth * ratio, phi_margin);
   [Gm, Pm, Wgm, Wpm] = margin(P * G * later.G);
   phi_diff = Pm - parameter.phimarginmin;
   phi_reg = min(phi_diff, phi_reg);   
   [~, phi_marginreg] = bode(P * G * later.G, bandwidth);
   phi_marginreg = min(45, 175 - abs(phi_marginreg));
   [trap, fval, exitflag] = trapdesign(P * later.G, G, bandwidth, num, phi_reg, Wpm, phi_creg, mag_creg, phi_marginreg);
   if exitflag == 1 || exitflag == 2
       bfailure = 1;
       break;
       %% ���Ч�����ҵ����õĽ�  ���ܵ��ţ� ��������ѭ��
        %% ά�����ڵĲ����� �޸�Լ��ʹ�����ܸ���
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

