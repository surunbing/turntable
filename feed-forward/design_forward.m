function [forward, exitflag] = design_forward(P, G, data, K, taum, taue, gain_max, bandwidth, maglim, philim, option)
%   给出前馈环节的增益

forward.G = 1;
forward.K = 0;
exitflag = 1;
kmax = gain_max;
kmin = 0;
k_pre = 0.5 * (kmax + kmin);
if strcmp(option.type, 'transfer-function')
    GC = (P * G + 1 / G * G * gain_max) / (1 + P * G);
    data_check = CL_check(GC, bandwidth, maglim, philim);
    if data_check.bmag ~= 1 || data_check.bphi ~= 1
        exitflag = 0;
    else
        %% 开始设计 二分法
        while 1
            k = 0.5 * (kmax + kmin);
            GC = (P * G + 1 / G * G * k) / (1 + P * G);
            data_check = CL_check(GC, bandwidth, maglim, philim);
            if data_check.bmag ~= 1 || data_check.bphi ~= 1
                kmin = k;
                kmax = kmax;
            else
                kmin = kmin;
                kmax = k;
                if abs((kmin + kmax) * 0.5 - k_pre) < 1e-3
                    k = (kmin + kmax) * 0.5;
                    break;
                end
                k_pre = k;
            end
        end
    end
elseif strcmp(option.type, 'discrete')
    
end
forward.G = 1 / G * k;
forward.K = k;
exitflag = 1;
end


function data_check = CL_check(GC, bandwidth, maglim, philim)
ncount = round(bandwidth / 2 / pi);
frequence = linspace(1, ncount, ncount) * 2 * pi;
data.mag = zeros(ncount, 1);
data.phi = zeros(ncount, 1);
data.fre = frequence;
[mag, phi] = bode(GC, data.fre);
for i = 1 : ncount
    data.mag(i) = mag(1, 1, i);
    data.phi(i) = phi(1, 1, i);
end
option.type = 'close-loop';
data_check = CLIndic_check(data, bandwidth, maglim, philim, option);
end

