close all;

parameter_init;

global parameter bflag_tradition;

%% ʹ��ʲô���
%% ���bandwidth
% if parameter.bandwidth >= 7 * 2 * pi &&  parameter.ratio <= 5
     bflag_tradition = 0;
     direct_test;
% else
%    bflag_tradition = 1;
 %   tradition_design;
% end
%     
