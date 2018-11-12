function [Rp] = GetRp(x, T)
%   »ñµÃĞ³Õñ·å
xi = x(1);
omegan = x(2);
if xi >= 1 / sqrt(2)
    Rp = 1;
else
% y = 1 - 2 * xi * xi;
% ar = 3 * T * T;
% br = 2 - 4 * T * T * omegan * omegan * (1 - 2 * xi * xi);
% deltar = (4 * T * T * omegan * omegan * y + 1) ^ 2 - 3 * (4 * T ^ 4 * omegan ^ 4 - 1);
% x1 = (-br + sqrt(deltar)) / (2 * ar);
% f = (omegan * omegan - x1) ^ 2 + 4 * xi * xi * omegan * omegan * x1;
% delta = (2 * T * T * omegan * omegan * (1 - 2 * xi * xi)) / 3 + 2 / 3;
% Rp = omegan * omegan / sqrt(f * (sqrt(deltar) / 6 + delta));
omega = sqrt(omegan * omegan - 2 * xi * xi * omegan * omegan);
Rp = omegan * omegan / sqrt((omegan * omegan - omegan * omegan) ^ 2 + (2 * xi * omegan *omega) ^ 2);
end
end

