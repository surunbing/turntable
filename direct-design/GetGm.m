function [kg] = GetGm(x, T)
%   获得幅值裕度
xi = x(1);
omegan = x(2);
kg = (2 * T * xi * omegan + 1) * (T * omegan + 2 * xi) / (T * omegan);
end

