function [LowGain] = GetPmloss(LowGain, wc)
%   ��ȡ����Ƶ�ʴ�����λ��ʧ
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

