clc, clear
close all

%% 补充低频增益，同时不破坏低频相位
%% 均使用弧度制, 角度制度

%% 定义结构体
LowGain.count = 0;
LowGain.fre = 0;
LowGain.alpha = 0;
LowGain.K = 0;
LowGain.Sum = 0;
LowGain.Pmloss = 0;

%% 选择频点，必须低频，尽量减少相位衰减对相角裕度的影响
%% 频点不能太低，否则影响真正的静态性能
LOW_FRE_MIN = 0.01 * 2 * pi;
LOW_FRE_MAX = 0.5 * 2 * pi;
LOW_PMLOSS_MAX = 3;

%% 如果低频增益影响相位裕度过大
%% 调整频带, 减小最高频带限制

