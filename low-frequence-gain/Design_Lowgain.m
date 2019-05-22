% load('control.mat');
%% 补充低频增益，同时不破坏低频相位
%% 均使用弧度制, 角度制度

%% 定义结构体
LowGain.count = 0;
LowGain.fre = 0;
LowGain.alpha = 0;
LowGain.K = 1;
LowGain.Pmloss = 0;

%% 选择频点，必须低频，尽量减少相位衰减对相角裕度的影响
%% 频点不能太低，否则影响真正的静态性能
LOW_FRE_MIN = 0.1 * 2 * pi;
LOW_FRE_MAX = 0.3 * 2 * pi;
LOW_PMLOSS_MAX = 1;

%% 如果低频增益影响相位裕度过大
%% 调整频带, 减小最高频带限制

LowGain = GetNewgain(LowGain, LOW_FRE_MIN, LOW_FRE_MAX, 3, 250);
% LowGain = GetNewgain(LowGain, LOW_FRE_MIN, LOW_FRE_MAX, 3, 250);
% LowGain = GetNewgain(LowGain, LOW_FRE_MIN, LOW_FRE_MAX, 3, 250);
% LowGain = GetNewgain(LowGain, LOW_FRE_MIN, LOW_FRE_MAX, 3, 250);
% LowGain = GetNewgain(LowGain, LOW_FRE_MIN, LOW_FRE_MAX, 3, 250);

Glow = GetlowgainG(LowGain);
% Glow = 1;
% LowGain.count = 0;
% figurename('迟后环节');
% bode(Glow);
% grid on

for i = 1 : LowGain.count
    tau = 1 / (sqrt(LowGain.alpha(i)) * LowGain.fre(i));
    LowGain.G(i) = tf([tau, 1], [LowGain.alpha(i) * tau, 1]);
end
