function [trap] = trapdesign(P, G, bandwidth)
%   ����ݲ�����
%   ������ ��ǰ����  ����ջ���Ƶ�������ջ���Ƶ�����������ҳ���һ�������ϵ�Ƶ�㣬���ͬ����ľ��룬
%   �Ƿ��ܶ����Ż����⣬Ѱ������Ϊ�ݲ��˲����������������޶��ڽ�Ϊխ�ķ�Χ�� must have a try
%   �������������ݲ��˲�������ô����
%   ָ�꺯��������

%% ���Ƶ��
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

%% ����һ�����ڣ� �Ż�ʵ��
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
% �õ��ٺ󻷽ڼ��� 
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

