function [wc] = GetWc(x, T)
%   »ñµÃ¼ôÇÐÆµÂÊ
xi = x(1);
omegan = x(2);
a = T * T;
b = 4 * xi * xi * T * T * omegan * omegan - 2 * T * T * omegan * omegan + 1;
c = omegan * omegan * (2 * xi + T *omegan) ^ 2;
d = - omegan ^ 4;
A = b * b - 3 * a * c;
B = b * c - 9 * a * d;
C = c * c - 3 * b * d;
M = (2 * A * b - 3 * a * B) / (2 * sqrt(A) ^ 3);
theta = acos(M);

delta = B * B - 4 * A * C;

if delta > 0 && A > 0
    xc = (-b - 2 * sqrt(A) * cos(theta / 3)) / (3 * a);
    wc = sqrt(xc);
elseif delta > 0 && A <= 0
    xc = (-b + sqrt(A) * (cos(theta / 3) - sqrt(3) * sin(theta / 3))) / (3 * a);
    wc = sqrt(xc);
elseif delta <= 0
    xc = (-b + sqrt(A) * (cos(theta / 3) + sqrt(3) * sin(theta / 3))) / (3 * a);
    wc = sqrt(xc); 
end
wc = abs(wc);
end

