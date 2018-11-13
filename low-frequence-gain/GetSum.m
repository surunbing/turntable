function [LowGain] = GetSum(LowGain)
%   获得低频总增益
LowGain.Sum = 1;
for i = 1 : count
   LowGain.Sum = LowGain.Sum * LowGain.K(i); 
end
end

