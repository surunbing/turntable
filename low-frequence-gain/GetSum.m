function [LowGain] = GetSum(LowGain)
%   ��õ�Ƶ������
LowGain.Sum = 1;
for i = 1 : count
   LowGain.Sum = LowGain.Sum * LowGain.K(i); 
end
end

