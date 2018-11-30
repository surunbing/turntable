function [trap, fval, exitflag] = trapdesign(P, G, bandwidth, num, pm_cost, wc, phi_reg, mag_reg)
%   ����ݲ�����
%   ������ ��ǰ����  ����ջ���Ƶ�������ջ���Ƶ�����������ҳ���һ�������ϵ�Ƶ�㣬���ͬ����ľ��룬
%   �Ƿ��ܶ����Ż����⣬Ѱ������Ϊ�ݲ��˲����������������޶��ڽ�Ϊխ�ķ�Χ�� must have a try
%   �������������ݲ��˲�������ô����
%   ָ�꺯��������
global parameter 
%% ���ƣ������Ȳ��˲��������￪ʼ������ָ��

%% ����һ�����ڣ� �Ż�ʵ��
ncount = round(bandwidth / 2 / pi);
data.fre = linspace(1, ncount, ncount)' * 2 * pi;
data.mag = zeros(ncount, 1);
data.phi = zeros(ncount, 1);
[mag, phi] = bode(P * G, data.fre);
for i = 1 : ncount
    data.mag(i) = mag(1, 1, i);
    data.phi(i) = phi(1, 1, i);
end

%% ���start point
start = zeros(num * 3, 1);
first_fre = bandwidth - 5;
for i = 1 : num
   start(i * 3 - 2) = 1.1;
   start(i * 3 - 1)= 15;
   start(i * 3) = first_fre - (i - 1) * 2 * pi; 
end

%% ���Լ��
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

