function [P, G, para] = direct_design()

global parameter

%% 考虑处理omegan的特征频率.
T = parameter.T;
Tmax = parameter.Tmax;
wc_max = parameter.wc_max;
div = parameter.Tdiv;
bandwidth = parameter.bandwidth;
flag = -1;
while flag ~= 9
    option.type = 'wc';
    [P, G, para, flag] = Getomegaxi(bandwidth, 0, T, option);
    if flag == 1
        if para.wc < wc_max
            %% 检查是否符合要求
            if T > Tmax && para.pm <= parameter.pmmax
                flag = 9;
            elseif T > Tmax && para.pm > parameter.pmmax
                flag = -1;
                T = max(T / div, Tmax);
            elseif T <= Tmax
                flag = 9;
            end
        else
            if T > Tmax
                T = max(T / div, Tmax);
            else
                bandwidth = bandwidth - pi;
            end
        end
    elseif flag ~= 1
        if T > Tmax
            T = max(T / div, Tmax);
        elseif T <= Tmax
            bandwidth = bandwidth - pi;
        end
    end
end
% P = P * tf(1, [1 / (100 * 2 * pi), 1]);
end

