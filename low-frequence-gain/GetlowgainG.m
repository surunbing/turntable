function [Glow] = GetlowgainG(LowGain)
%   获得低频增益传递函数
Glow = 1;
for i = 1 : LowGain.count
    tau = 1 / (sqrt(LowGain.alpha(i)) * LowGain.fre(i));
    Glow = Glow * tf([tau, 1], [LowGain.alpha(i) * tau, 1]);
end
Glow = Glow * LowGain.K;
end

