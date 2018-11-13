function [LowGain] = GetPmloss(LowGain, wc)
%   获取剪切频率处的相位损失
LowGain.Pmloss = 0;
for i = 1 : count
    alpha = LowGain.alpha(i);
    fre = LowGain.fre(i);
    tau = 1 / (sqrt(alpha) * fre);
    T = alpha * tau;
    G = tf([tau, 1], [T, 1]);
    [~, phi] = bode(G, wc);
    LowGain.Pmloss = LowGain.Pmloss + phi;
end
end

