function [c, ceq] = nlconws(x, T, pmr, kgr, Mre, wcmax)
%UNTITLED10 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
c(1) = pmr - GetPm(x, T);
c(2) = kgr - 20 * log10(GetGm(x, T));
c(3) = 20 * log10(GetRp(x, T)) - Mre;
c(4) = GetWc(x, T) - wcmax;
ceq = 0;
end

