function [kg] = GetGm(x, T)
%   ��÷�ֵԣ��
xi = x(1);
omegan = x(2);
kg = (2 * T * xi * omegan + 1) * (T * omegan + 2 * xi) / (T * omegan);
end

