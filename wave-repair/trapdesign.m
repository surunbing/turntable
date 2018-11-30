function [trap, fval, exitflag] = trapdesign(P, G, bandwidth, num, pm_cost, wc, phi_reg, mag_reg)
%   设计陷波环节
%   如何设计 从前向后加  假设闭环相频单减，闭环幅频单调递增，找出第一个不符合的频点，检查同带宽的距离，
%   是否能定义优化问题，寻优内容为陷波滤波器的三个参数，限定在较为窄的范围内 must have a try
%   如果是添加两个陷波滤波器，怎么考虑
%   指标函数如何添加
global parameter 
%% 改善，不加先波滤波器，哪里开始不满足指标

%% 加入一个环节， 优化实验
ncount = round(bandwidth / 2 / pi);
data.fre = linspace(1, ncount, ncount)' * 2 * pi;
data.mag = zeros(ncount, 1);
data.phi = zeros(ncount, 1);
[mag, phi] = bode(P * G, data.fre);
for i = 1 : ncount
    data.mag(i) = mag(1, 1, i);
    data.phi(i) = phi(1, 1, i);
end

%% 获得start point
start = zeros(num * 3, 1);
first_fre = bandwidth - 5;
for i = 1 : num
   start(i * 3 - 2) = 1.1;
   start(i * 3 - 1)= 15;
   start(i * 3) = first_fre - (i - 1) * 2 * pi; 
end

%% 获得约束
lb = zeros(num * 3, 1);
ub = zeros(num * 3, 1);
for i = 1 : num
    lb(i * 3 - 2) = 1;
    lb(i * 3 - 1) = parameter.trapTmin;
    lb(i * 3) = parameter.trapfremin;
    ub(i * 3 - 2) = 5;
    ub(i * 3 - 1)= 20;
    ub(i * 3) = parameter.trapfremax;
end

% options = optimset('Algorithm','interior-point');
options = optimset('Algorithm','sqp');

[x, fval, exitflag] = fmincon(@(x)GetTrapcost(x, num, data, 0, 0)...
    , start, [], [], [], [], lb, ub, @(x)nonlcon_trap(x, wc, pm_cost, data, num, phi_reg, mag_reg), options);
toc
% [c, ceq] = nonlcon_trap(x, wc, pm_cost, data);
trap.num = num;
trap.e = zeros(num, 1);
trap.T = zeros(num, 1);
trap.f = zeros(num, 1);
% trap.G = 0;
for i = 1 : num
    trap.e(i) = x(i * 3 - 2);
    trap.T(i) = x(i * 3 - 1);
    trap.f(i) = x(i * 3);
    trap.G(i) = tf([1, trap.e(i) * trap.T(i), trap.f(i) * trap.f(i)], [1, trap.T(i), trap.f(i) * trap.f(i)]);
end

end

