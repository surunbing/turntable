function [bStable, bGm, bPm, bPhi, bWc] = Stability_check(data, data2, Gmmin, Pmmin, phi_min, wc_max, bandwidth, ratio, option)
%   检查系统的稳定性， 判断幅值裕度和相角裕度是否满足要求，剪切频率 判断是否条件稳定 判断
if strcmp(option.type, 'transfer') == 1
    bStable = isstable(data);
    [Gm, Pm, ~, Wc] = margin(data);
    frequence = logspace(log10(0.1), log10(bandwidth * ratio), 100);
    [mag, phi] = bode(data, frequence);
    Mag = zeros(length(frequence), 1);
    Phi = zeros(length(frequence), 1);
    for i = 1 : length(frequence)
       Mag(i) = mag(1, 1, i);
       Phi(i) = phi(1, 1, i);
    end
    phimin = min(Phi);
    if phimin < phi_min
        bPhi = 0;
    else
        bPhi = 1;
    end
    
    Gm = 20 * log10(Gm);
    if Gm < Gmmin
        bGm = 0;
    else
        bGm = 1;
    end
    
    if Pm > Pmmin
        bPm = 1;
    else
        bPm = 0;
    end
    if Wc > wc_max
        bWc = 0;
    else
        bWc = 1;
    end
elseif strcmp(option.type, 'discrete') == 1
    %% 判断穿越频率
    ncount = data.fre;
    for i = 1 : ncount
        
    end
end
end

