function [flag, advance, data_cur] = DesignPhi(data)
%DESIGNPHI 此处显示有关此函数的摘要
%   此处显示详细说明
global START_RATIO PHI_ADVANCE PHI_ADVANCE_GARGIN BANDWIDTH DELTA_RATIO RATIO_MIN
ratio = START_RATIO;
while 1
    advance = FillPhase(data, ratio, PHI_ADVANCE, PHI_ADVANCE_GARGIN, BANDWIDTH);
    bsucceed = checkFillPhase(data, advance);
    if bsucceed == 1
        break;
    else
       %%缩小剪切频率
       if ratio < RATIO_MIN
           flag = 0;
       else
           ratio = ratio - DELTA_RATIO;
           if ratio < 2.5
               ratio = 2;
           end
       end
    end          
end
frequence = data.fre;
[mag_P, phi_P] = bode(advance.P, frequence);
mag = zeros(length(frequence), 1);
phi = zeros(length(frequence), 1);
for i = 1 : length(frequence)
    mag(i) = 20 * log10(mag_P(1, 1, i));
    phi(i) = phi_P(1, 1, i);
end
data_cur.fre = data.fre;
data_cur.mag = mag + data.mag;
data_cur.phi = phi + data.phi;
end

