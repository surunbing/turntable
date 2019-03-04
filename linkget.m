close all;

parameter_init;

global parameter bflag_tradition;

%% 使用什么设计
%% 检查bandwidth
% if parameter.bandwidth >= 7 * 2 * pi &&  parameter.ratio <= 5
%     bflag_tradition = 0;
%     direct_test;
% else
    bflag_tradition = 1;
    tradition_design;
% end
%     
% e = 0.92;
% T = 20;
% f = 10 * 2 * pi;
% trap = tf([1 e * T f * f], [1, T f * f]);
% bode(trap);
% grid on