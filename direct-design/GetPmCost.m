function [cost] = GetPmCost(x, T)
%   获得代价函数  -pm
cost = -GetPm(x, T);
end

