function [cost] = GetWsCost(x, T)
%   获得代价函数  -pm
cost = -GetDt(x, T);
end

