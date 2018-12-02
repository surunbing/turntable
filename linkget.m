close all;

parameter_init;

global parameter

%% ¼ì²ébandwidth
if parameter.bandwidth >= 7 * 2 * pi &&  parameter.ratio <= 5
    direct_test;
else
    tradition_design;
end
    