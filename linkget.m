close all;

parameter_init;

global parameter bflag_tradition;

%% 使用什么设计
%% 检查bandwidth
if parameter.bandwidth >= 3 * 2 * pi &&  parameter.ratio >= 2 && methodmodel == 0
    bflag_tradition = 1;
    tradition_design;
elseif parameter.bandwidth >= 7 * 2 * pi &&  parameter.ratio <= 2 && methodmodel == 0
    bflag_tradition = 0;
    direct_test;
elseif methodmodel == 1
    bflag_tradition = 0;
    direct_test;
elseif methodmodel == 2
    bflag_tradition = 1;
    tradition_design;
end
