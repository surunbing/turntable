function [advance, ratio] = FillPhase(data, ratio, phi_advance, phi_advance_margin, bandwidth)

%补相位
advance.phi_advance = phi_advance;
advance.phi_advance_margin = phi_advance_margin;
% 穿越频率3-5倍
fre_through = bandwidth * ratio;
if fre_through > 380
    fre_through = 380;
    ratio = fre_through / bandwidth;
end
fre_error = abs(data.fre - fre_through * ones(1, length(data.fre)));
fre = min(fre_error);
num = find(fre_error == fre);
advance.fre = data.fre(num);
advance.phi = data.phi(num);
advance.phi_error = abs(advance.phi) - 135; 
advance.num = num;

%%拆分超前
num1 = advance.phi_error / phi_advance;
advance.count = fix(num1);
advance.phi_remain = advance.phi_error - phi_advance * advance.count;
if(advance.phi_remain < 0.5)
    advance.phi_remain = 0;
end

%% 设计环节
m = sin((phi_advance + phi_advance_margin) / 180 * pi);
alpha = (1 - m) / (1 + m);
tau = 1 / (sqrt(alpha) * advance.fre);
T = alpha * tau;
advance.G_advance = tf([tau, 1], [T, 1]);

if(advance.phi_remain ~= 0)
    m1 = sin((advance.phi_remain + phi_advance_margin) / 180 * pi);
    alpha1 = (1 - m1) / (1 + m1);
    tau1 = 1 / (sqrt(alpha1) * advance.fre);
    T1 = alpha1 * tau1;
    advance.G_advance1 = tf([tau1, 1], [T1, 1]);
else
    advance.G_advance1 = 1;
end

P = 1;
for i = 1: advance.count
    P = P * advance.G_advance;
end
P = P * advance.G_advance1;
%% 求取增益
[mag, ~] = bode(P, data.fre(num));
Mag = 20 * log10(mag(1, 1, 1)) + data.mag(num);
advance.gain = 1 / 10 ^ (Mag / 20);    
advance.P = P;

end

