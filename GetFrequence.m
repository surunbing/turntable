function [fre_Rp, fre_con] = GetFrequence(bandwidth, ratio, nRp, nCon, nCRp)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明

fre_start = log10(bandwidth * ratio * 0.8);
fre_end = log10(bandwidth * ratio * 1.7);
fre_Rp = logspace(fre_start, fre_end,  nRp)';

fre_start = log10(bandwidth * 0.3);
fre_end = log10(bandwidth * ratio * 0.9);
fre_con = [logspace(fre_start, fre_end,  nCon), bandwidth * 3, bandwidth * 5, logspace(2, 3,  nCRp)]';

end

