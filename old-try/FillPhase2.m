function [advance] = FillPhase(data, ratio, phi_advance, phi_advance_margin, bandwidth)
global parameter
%补相位

advance.phi_advance = phi_advance;
advance.phi_advance_margin = phi_advance_margin;
% 穿越频率3-5倍
fre_through = bandwidth * ratio;
if fre_through > 380
    fre_through = 380;
    ratio = fre_through / bandwidth;
end
fre_error = abs(data.fre - fre_through * ones(length(data.fre), 1));
fre = min(fre_error);
num = find(fre_error == fre);
advance.fre = data.fre(num);
advance.phi = data.phi(num);
advance.phi_error = abs(advance.phi) - parameter.tradphi; 
advance.num = num;

start_n = 2;
if advance.phi_error > 800 && advance.phi_error < 180
    start_n = 2;
elseif advance.phi_error <= 80
    start_n = 1;
else
    start_n = 3;
end

cost = zeros(10 - start_n + 1, 1);

for n = start_n : 10
    angle_margin = advance.phi_error / 180 * pi;
    a  = (1 - sin(angle_margin / n)) / (1 + sin(angle_margin / n));
    wm = advance.fre;
    w = wm * 3;
    ccost = sqrt(w * w + a * wm * wm) / sqrt(a * a * w * w + a * wm * wm);
    ccost = 20 * log10(ccost);
    cost(n - start_n + 1) = ccost + n;
end
cost_min = min(cost);
num_min = find(cost == cost_min);


%%拆分超前
num1 = num_min(1) + start_n - 1;
advance.count = fix(num1);
advance.phi_remain = 0;
if(advance.phi_remain < 0.5)
    advance.phi_remain = 0;
end

%% 设计环节
m = sin((advance.phi_error / num1) / 180 * pi);
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
    advance.G(i) = advance.G_advance;
end
advance.G(advance.count + 1) = advance.G_advance1;
P = P * advance.G_advance1;
%% 求取增益
[mag, ~] = bode(P, data.fre(num));
Mag = 20 * log10(mag(1, 1, 1)) + data.mag(num);
advance.gain = 1 / 10 ^ (Mag / 20);    
advance.P = P * advance.gain;
advance.ratio = ratio;

end

