function [trap] = trapdesign(P, G, bandwidth)
%   设计陷波环节
%   如何设计 从前向后加  假设闭环相频单减，闭环幅频单调递增，找出第一个不符合的频点，检查同带宽的距离，
%   是否能定义优化问题，寻优内容为陷波滤波器的三个参数，限定在较为窄的范围内 must have a try
%   如果是添加两个陷波滤波器，怎么考虑
%   指标函数如何添加

%% 检查频点
% ncount = round(bandwidth / 2 / pi);
% data.mag = zeros(ncount, 1);
% data.phi = zeros(ncount, 1);
% data.fre = linspace(1, ncount) * 2 * pi;
% [mag, phi] = bode(P * G / (1 + P * G), data.fre);
% for i = 1 : 1 : ncount
%    data.mag(i) = mag(1, 1, i);
%    data.phi(i) = phi(1, 1, i);
% end
% option.type = 'close-loop';
% data_check = CLIndic_check(data, bandwidth, option);

%% 加入一个环节， 优化实验
[~, ~, ~, wc] = margin(P * G);
pm_cost = 10;
ncount = round(bandwidth / 2 / pi);
data.fre = linspace(1, ncount, ncount) * 2 * pi;
data.mag = zeros(ncount, 1);
data.phi = zeros(ncount, 1);
[mag, phi] = bode(P * G, data.fre);
for i = 1 : ncount
    data.mag(i) = mag(1, 1, i);
    data.phi(i) = phi(1, 1, i);
end

tic
% 得到迟后环节计算 
start = [1.1, 5, 100, 1.1, 5, 90];
lb = [1; 10; bandwidth * 0.5; 1; 4; bandwidth * 0.5];
ub = [5; 70; bandwidth + 10; 5; 70; bandwidth + 10];
% options = optimset('Algorithm','interior-point');
options = optimset('Algorithm','sqp');
[x, fval, exitflag] = fmincon(@(x)GetTrapcost(x, 0, data, 0, 0)...
    , start, [], [], [], [], lb, ub, @(x)nonlcon_trap(x, wc, pm_cost, data), options);
toc
[c, ceq] = nonlcon_trap(x, wc, pm_cost, data);
trap.G1 = tf([1, x(1) * x(2), x(3) * x(3)], [1, x(2), x(3) * x(3)]);
trap.G2 = tf([1, x(4) * x(5), x(6) * x(6)], [1, x(5), x(6) * x(6)]);
end

