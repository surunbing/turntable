function [Design] = Designplastic(data, advance, data_cur)
%   波形整形 使用超前迟后环节和非典型的二阶环节   形状判定方法可以改进
%   约束最小相位， 计算cost function
global PHI_MIN BANDWIDTH

%% 检查是否可以进行设计
%% 肯定可以进行设计
num = find(data.fre == advance.fre);
frequence = data.fre(1 : num);
phi_min = min(data_cur.phi);
num = find(data_cur.phi == phi_min);
num = num(1);
fre_phimin = frequence(num);
if fre_phimin < -90
    %% 在频率后最大相位与最大频点
    num = find(data_cur.fre > 2 * pi);
    num = num(1);
    phi_max = data_cur.phi(num : end);
    num = find(data_cur.phi == phi_max);
    num = num(1);
    fre_phimax = data_cur.fre(num);
    if(fre_phimax < BANDWIDTH)
        bPhiFillSlow = 1;
    else 
        bPhiFillSlow = 0;
    end
    %% 下缓慢则标记，并且记录
    if bPhiFillSlow == 1
        num = find(data_cur.fre > BANDWIDTH);
        num = num(1);
        phi_max = data_cur.phi(num : end);
        num = find(data_cur.phi == phi_max);
        num = num(1);
        fre_phimax = data_cur.fre(num);
    end
    phi_reg = fre_phimin - PHI_MIN;
    % 正常整形;
    bNormal = 1;
else
    % 不可以正常设计
    bNormal = 0;
end




% 获得闭环
end

