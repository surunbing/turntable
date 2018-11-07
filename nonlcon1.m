% G ����
function [c,ceq] = nonlcon1(x, data, series, P, bandwidth, nCon, nCRp, phi_reg, mag_max)

[mag, phi] = GetMagPhi(x, series, data.fre);
mag = mag + data.mag;
phi = phi +data.phi;
nCount = nCon;
c = zeros(nCon + 2 + nCRp, 1);
for i = 1 : nCon
    c(i) = - phi(i) - phi_reg;
end
c(nCount + 1) = -mag(nCount + 1);
c(nCount + 2) = mag(nCount + 2);
% �ջ�����
Mag = mag;
Phi = phi;
complex_bode = 10 .^ (Mag ./ 20) .* complex(cos(Phi ./ 180 .* pi), sin(Phi ./ 180 .* pi));
% �ջ�
complex_bode = complex_bode ./ (1 + complex_bode);
%ת��bode��ʽ
Mag = log10(abs(complex_bode)) * 20;
%Mag_max = max(Mag);
for i = 1 : nCRp
    c(nCount + 2 + i) = Mag(nCon + 2 + i) - mag_max;
end
ceq = [];
