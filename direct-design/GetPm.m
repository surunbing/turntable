function [Pm] = GetPm(x,T)
%   获得相位裕度
xi = x(1);
omegan = x(2);
wc = GetWc(x, T);
v1 = wc * (2 * T * xi * omegan + 1) / (omegan * (T * omegan + 2 * xi) - T * wc * wc);
v2 = atan(v1) / pi * 180;
Pm = 90 - v2;
end

