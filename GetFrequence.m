function [fre_Rp, fre_con] = GetFrequence(bandwidth, ratio)
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

fre_start = log10(bandwidth * ratio * 0.8);
fre_end = log10(bandwidth * ratio * 1.7);
fre_Rp = logspace(fre_start, fre_end,  5)';

fre_start = log10(bandwidth);
fre_end = log10(bandwidth * ratio * 0.9);
fre_con = [logspace(fre_start, fre_end,  5), bandwidth * 3, bandwidth * 5]';

end

