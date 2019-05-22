function [trap] = design_instruction_preprocessing(data)
global parameter

phi_min = min(data.phi);
phi_min = phi_min(1);
maglim = parameter.maglim;
if phi_min > -10 && phi_min < -parameter.philim 
    philim = (10 + abs(phi_min)) / 2;
else
    philim = parameter.philim;
end

%% 加入一个环节， 优化实验
bandwidth = parameter.bandwidth;


for num = 1 : 3
    %% 获得start point
    start = zeros(num * 3, 1);
    first_fre = bandwidth - 5;
    for i = 1 : num
       start(i * 3 - 2) = 1;
       start(i * 3 - 1)= 15;
       start(i * 3) = 50; 
    end

    %% 获得约束
    lb = zeros(num * 3, 1);
    ub = zeros(num * 3, 1);
    for i = 1 : num
        lb(i * 3 - 2) = 0.1;
<<<<<<< HEAD
        lb(i * 3 - 1) = parameter.trapTmin;
        lb(i * 3) = parameter.trapfremin;
        ub(i * 3 - 2) = 1.01;
        ub(i * 3 - 1)= 500;
        ub(i * 3) = parameter.trapfremax + 200 ;
=======
        lb(i * 3 - 1) = 0.1 * parameter.trapTmin;
        lb(i * 3) = parameter.trapfremin;
        ub(i * 3 - 2) = 5;
        ub(i * 3 - 1)= 500;
        ub(i * 3) = parameter.trapfremax * 2;
>>>>>>> 26649d5039a8c3eaddce917adfb24bbcfa761eb2
    end

    %  options = optimset('Algorithm','interior-point');
    options = optimset('Algorithm','sqp', 'MaxFunEvals', 2500);

    [x, fval, exitflag] = fmincon(@(x)GetTrapcostpre(x, num, data)...
        , start, [], [], [], [], lb, ub, @(x)nonlcon_trappre(x, data, num, philim, maglim), options);
    toc
    if exitflag == 1 || exitflag == 2
        break;
    end
end
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

