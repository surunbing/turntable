function [q] = myfun(x, xi, omegan, T)
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
para1 = omegan * omegan;
para2 = complex(0, x);
para3 = complex(omegan * (T * omegan + 2 * xi) - T * T * x * x, 2 * (T * xi * omegan + 1) * x);
complex_wc = para1 / (para2 * para3);
q = angle(complex_wc) / pi * 180 + 130;
end

