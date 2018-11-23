function [trap, later, bfailure, data_check] = wave_repair(P, G, para, wc_up, data, bandwidth)
%   repair ���Σ�ͨ������wc�� �����ǵ�ԣ�ȣ� trap���� ��Լ����Χ���������ε�˫ʮָ����
ratio_max = wc_up / bandwidth;
%% �Ƿ���Ҫ�����ܷ���Ƴ�������
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
       %% ���Ч�����ҵ����õĽ�
        %% ά�����ڵĲ����� �޸�Լ��ʹ�����ܸ���
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



% figurename('�ٺ�');
% margin(P * G * later.G);
% grid on
% figurename('�ٺ�ջ�');
% bode(P * G * later.G / (1 + P * G * later.G));
% grid on

% num = 4;
% K = P * G * later.G;
% for i = 1 : num
%    K = K * trap.G(i); 
% end
% figurename('�ݲ��˲���');
% margin(K);
% grid on
% figurename('�ݲ��˲����ջ�');
% bode(K / (1 + K));
% grid on
end

