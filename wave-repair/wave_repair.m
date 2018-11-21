function [trap, later, flag, data_check] = wave_repair(P, G, para, wc_up, data, bandwidth)
%   repair 波形，通过调整wc， 期望角的裕度， trap个数 与约束范围来调整波形到双十指标内


[later, fval, exitflag] = Holddonewc(P, G, para, data, wc_r, phi_dist, option)


outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

