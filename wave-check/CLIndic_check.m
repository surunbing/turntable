function [data_check] = CLIndic_check(data, bandwidth, option)
% 检查双十特性 开环检查或闭环检查
ncount = round(bandwidth / 2 / pi);
data_check.mag = zeros(ncount, 1);
data_check.phi = zeros(ncount, 1);
data_check.mag11_reg = zeros(ncount, 1);
data_check.mag09_reg = zeros(ncount, 1);
data_check.phips_reg = zeros(ncount, 1);
data_check.phing_reg = zeros(ncount, 1);
data_check.bmag = 1;
data_check.bphi = 1;
if strcmp(option.type, 'close-loop') == 1
    for i = 1 : ncount
       if data.mag(i) > 0.9 && data.mag(i) < 1.1
           data_check.mag(i) = 1;
       else
           data_check.mag(i) = 0;
           data_check.bmag = 0;
       end
       if abs(data.phi(i)) < 8
           data_check.phi(i) = 1;
       else
           data_check.phi(i) = 0;
           data_check.bphi = 0;
       end
       %% 计算差值
       data_check.mag11_reg(i) = 1.1 - data.mag(i);
       data_check.mag09_reg(i) = data.mag(i) - 0.9;
       data_check.phips_reg(i) = 10 - data.phi(i);
       data_check.phing_reg(i) = data.phi(i) + 10;
    end
elseif strcmp(option.type, 'open-loop') == 1
    %% 由开环计算闭环
    data_close = data;
    close_complex = data.mag .* complex(cos(data.phi / 180 * pi), sin(data.phi / 180 * pi));
    close_complex = close_complex ./ (1 + close_complex);
    data_close.mag = abs(close_complex);
    data_close.phi = angle(close_complex) / pi * 180;
    for i = 1 : ncount
        if data_close.mag(i) > 0.9 && data_close.mag(i) < 1.1
           data_check.mag(i) = 1;
       else
           data_check.mag(i) = 0;
           data_check.bmag = 0;
       end
       if abs(data_close.phi(i)) < 9
           data_check.phi(i) = 1;
       else
           data_check.phi(i) = 0;
           data_check.bphi = 0;
       end
       %% 计算差值
       data_check.mag11_reg(i) = 1.1 - data_close.mag(i);
       data_check.mag09_reg(i) = data_close.mag(i) - 0.9;
       data_check.phips_reg(i) = 10 - data_close.phi(i);
       data_check.phing_reg(i) = data_close.phi(i) + 10;
    end
end
end

