% G ����
function [c,ceq] = nonlcon1(x, data, series, P, bandwidth, nCon, nCRp)

[mag, phi] = GetMagPhi(x, series, data.fre);
mag = mag + data.mag;
phi = phi +data.phi;
nCount = nCon;
c = zeros(nCon + 4 + nCRp, 1);
for i = 1 : nCon
    c(i) = - phi(i) - 175;
end
c(nCount + 1) = -mag(nCount + 1);
c(nCount + 2) = mag(nCount + 2);
c(nCount + 3) =  - phi(nCount + 1) - 180;
c(nCount + 4) =  - phi(nCount + 2) - 190;
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
    c(nCount + 4 + i) = Mag(nCon + 2 + i) - 6;
end
ceq = [];
