function [flag, advance] = DesignPhi(data, bandwidth)
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
           
end

end

