function [later] = Boostedlfgain(data, advance, fre, alpha)
%UNTITLED3����Ƶ����
%   fre: ����Ƶ��
[mag, ~] = bode(advance.P * advance.gain, advance.fre);
Mag_GP = 20 * log10(mag) + data.mag(advance.num);      %% ��������λԣ�ȳ�������
tau = 1 / (sqrt(alpha) * fre);
T = alpha * tau;
later.G_later = tf([tau, 1], [T, 1]);
later.fre = fre;
later.alpha = alpha;
[mag, ~] = bode(later.G_later, advance.fre);
mag = 20 * log10(mag) + Mag_GP;
later.gain = 1 / (10 ^ (mag / 20));
end

