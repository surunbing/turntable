function [forward, exitflag] = design_forward(P, option, philim)
%   给出前馈环节的增益
global parameter

K = parameter.K;
taum = parameter.taum;
taue = parameter.taue;

forward.G = 1;
forward.K = 0;
exitflag = 1;
maglim = 2.0;%parameter.maglim;
% philim = parameter.philim;
kmax = parameter.forwardKmax;
bandwidth = parameter.bandwidth;
kmin = 0;
k_pre = 0.5 * (kmax + kmin);
G_auxiliary = tf(1, conv([1 / (parameter.para_aux1 * 2 * pi), 1], [1 / (parameter.para_aux2 * 2 * pi), 1]));

G = tf(K, [taue * taum, taue + taum, 1, 0]);

ts = parameter.Ts;
forward_n = 10;
G_speed = 1 / ts / forward_n / forward_n;
G_speed_sum = 0;
for i = 1 : forward_n
    G_speed_sum = G_speed_sum + tf(1, [(i - 1) * ts, 1]) - tf(1, [(i + 9) * ts, 1]);
end
G_speed = G_speed * G_speed_sum;

G_model = tf(K, [taue * taum, taue + taum, 1]);
G_forward = G_auxiliary / G_model * G_speed;
k = 0;
if strcmp(option.type, 'transfer-function')
    GC = (P * G + G_forward * G * parameter.forwardKmax) / (1 + P * G);
    data_check = CL_check(GC, bandwidth, maglim, philim);
    if data_check.bmag ~= 1 || data_check.bphi ~= 1
        exitflag = 0;
        k = parameter.forwardKmax;
    else
        %% 开始设计 二分法
        while 1
            k = 0.5 * (kmax + kmin);
            GC = (P * G + G_forward * G * k) / (1 + P * G);
            data_check = CL_check(GC, bandwidth, maglim, philim);
            if data_check.bmag ~= 1 || data_check.bphi ~= 1
                kmin = k;
%                 kmax = kmax;
            else
%                 kmin = kmin;
                kmax = k;
                if abs((kmin + kmax) * 0.5 - k_pre) < 1e-3
                    k = (kmin + kmax) * 0.5;
                    exitflag = 1;
                    break;
                end
                k_pre = k;
            end
        end
    end
elseif strcmp(option.type, 'discrete')
    
end
forward.G = G_forward * k;
forward.K = k;
% exitflag = 1;
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
    if data.phi(i) > 180
        data.phi(i) = data.phi(i) - 360;
    end
end
option.type = 'close-loop';
data_check = CLIndic_check(data, bandwidth, 10 ^ (maglim / 20) - 1, philim, option);
end

