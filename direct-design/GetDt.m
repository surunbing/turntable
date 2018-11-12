function [dt] = GetDt(x, T)
%   获得双十
xi = x(1);
omegan = x(2);

%% 相关变量
y = 1 - 2 * xi * xi;
as = T * T;
bs = 1 - 2 * y * as * omegan * omegan;
cs = omegan * omegan * (T * T * omegan * omegan - 2 * y);
ds = 21 * omegan ^ 4 / 121;
As = bs * bs - 3 * as * cs;
Bs = bs * cs - 9 * as * ds;
Cs = cs * cs - 3 * bs * ds;
Ts = (2  * As * bs - 3 * as * Bs) / (2 * sqrt(As) ^ 3);
thetas = acos(Ts);

xa11 = (-bs + sqrt(As * (cos(thetas / 3) - sqrt(3) * sin(thetas / 3)))) / (3 * as);

if xi > 1 / sqrt(2)
    omegaa11 = inf;
else
    omegaa11 = sqrt(xa11);
end

omegap10 = (sqrt((T * omegan + 2 * xi) ^ 2 + 8 * T * xi * omegan * tan(pi / 18) ^ 2) - (T * omegan + 2 * xi)) / (4 * T * xi * tan(pi / 18));

if omegaa11 > omegap10
    dt = omegap10;
elseif omegaa11 <= omegap10
    dt = omegaa11;
end

end

