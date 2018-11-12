function [bsucceed] = checkFillPhase(data, advance)
global PHI_MIN
%   判定是否可以整形: 补充相位是否可行， 
%% checkout if is Conditional stability
num = find(data.fre == advance.fre);
frequence = data.fre(1 : num);
[mag_P, phi_P] = bode(advance.P, frequence);
mag = zeros(length(frequence), 1);
phi = zeros(length(frequence), 1);
for i = 1 : length(frequence)
    mag(i) = 20 * log10(mag_P(1, 1, i));
    phi(i) = phi_P(1, 1, i);
end
% Mag = mag + data.mag(1 : length(frequence));
Phi = phi + data.phi(1 : length(frequence));
phi_min = min(Phi);
num = find(Phi == phi_min);
num = num(1);
bsucceed = 0;
fre_phimin = data.fre(num);
if fre_phimin < data.fre && fre_phimin > data.fre(1) && fre_phi_min > (PHI_MIN + 3)
    bsucceed = 1;
else
    bsuccees = 0;
end
end

