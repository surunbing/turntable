function [LowGain] = GetNewgain(LowGain, LOW_FRE_MIN, LOW_FRE_MAX, K, wc)
%% 调整所有频点，并且重新计算alpha和fre
count = LowGain.count + 1;
frequence = logspace(log10(LOW_FRE_MIN), log10(LOW_FRE_MAX), count);
LowGain.fre = frequence(1 : count);
LowGain.count = count;
LowGain.K = LowGain.K * K;
LowGain.alpha = zeros(count, 1) + LowGain.K ^ (1 / count);
LowGain = GetPmloss(LowGain, wc);
end

