function [LowGain] = GetNewgain(LowGain, LOW_FRE_MIN, LOW_FRE_MAX, K)
%% ��������Ƶ�㣬�������¼���alpha��fre
count = LowGain.count + 1;
frequence = logspace(log10(LOW_FRE_MIN), log10(LOW_FRE_MAX), count + 2);
LowGain.fre = frequence(2 : count + 1);
LowGain.count = count;
LowGain.K = LowGain.K * K;
LowGain.alpha = zeros(count, 1) + LowGain.K / count;
LowGain = GetSum(LowGain);
LowGain = GetPmloss(LowGain, wc);
end

