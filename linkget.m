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
close all
e = 0.90;
T = 4;
f = 0.5 * 2 * pi;
e1 = 0.90;
T1 = 4;
f1 = 0.1 * 2 * pi;
trap = tf([1 e * T f * f], [1, T f * f]);
trap1 = tf([1 e1 * T1 f1 * f1], [1, T1 f1 * f1]);
bode(trap);
hold on
bode(trap1);
hold on
bode(trap * trap1);
grid on