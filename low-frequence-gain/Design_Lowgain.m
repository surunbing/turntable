clc, clear
close all

%% 补充低频增益，同时不破坏低频相位

%% 定义结构体
LowGain.count = 0;
LowGain.fre = 0;
LowGain.alpha = 0;
LowGain.K = 0;
LowGain.Sum = 0;