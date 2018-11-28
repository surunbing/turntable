function [P, G, para] = direct_design(bandwidth, wc_max, K, taum, taue)
%% ���Ǵ���omegan������Ƶ��.
T = 0.0014 * 3;
Tmax = 0.0014 / 1.6;
flag = -1;
div = 1.05;
while flag ~= 9
    option.type = 'wc';
    [P, G, para, flag] = Getomegaxi(bandwidth, 0, T, K, taum, taue, option);
    if flag == 1
        if para.wc < wc_max
            %% ����Ƿ����Ҫ��
            if T > Tmax && para.pm <= 45
                flag = 9;
            elseif T > Tmax && para.pm > 45
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
end

