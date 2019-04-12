bandwidth = 10; % 带宽设置, hz 请设置整数倍
ratio = 3.0;

global parameter

%% 模型参数
% K = 1.096105440623064e+02;
% taue = 0.002143513958299;
% taum = 0.190856829212321;

% K = 496.7296;    %% ET205A
% taue = 0.0019;
% taum = 2.0624;

% K = 434;
% taum = 0.67;
% taue = 0.0035;

K = 1.56 * 180 / pi;
taue = 0.0039035;
taum = 0.984871194396488;

% nZLYCLStart = 0;
% nZLYCLEnd = 3;
% nQKStart = 4;
% nQKEnd = 6;
% nJZStart = 7;
% nJZEnd = 29;
% nDOF = 0;
% 
% outputfile = 'C:\Users\Momenta\Documents\毕业设计\助手\CSDA_FANGXUN\turntable/controller5.txt';

parameter.K = K;
parameter.taum = taum;
parameter.taue = taue;

parameter.bandwidth = bandwidth * 2 * pi;
parameter.ratio = ratio;        % 选择穿越频率倍数
parameter.wc_max = 650;     % 直接设计的时候最大带宽限制，请勿修改
parameter.T = 0.0014;       % 直接设计T
parameter.Tmax = parameter.T / 1.6; % 直接设计的最大T 
parameter.Tdiv = 1.05; %直接涉及惯性系数搜索因子
parameter.kgr = 5;     % 直接涉及相位裕度限制
parameter.Mre = 6;     % 直接设计谐振峰限制
parameter.pmr = 35;    % 直接设计相位裕度限制
parameter.dampmax = 1; % 直接设计阻尼比上限 
parameter.pmmax = 45;

%% 波形整型的相关参数
parameter.phi_creg = 8;   % 闭环整形相位目标
parameter.mag_creg = 0.8; % 闭环整形幅值目标, db
parameter.maglim = 0.7;  % 闭环最终幅值 db
parameter.philim = 8;     % 闭环最终相位
parameter.num_max = 3;    % 整型最大非线性环节数量
parameter.rdiv = 1.05;    % 指标优化或放开的系数
parameter.start_ratio = parameter.ratio * 0.8;
parameter.phimarginmin = 44.5;  % 相位裕度最小值

parameter.phi_margin = 122;  %闭环整形为了使用险波环节, 需要给出优化的相位裕度    这个数值可以寻优 但经过实验，效果不大
parameter.phi_reg = 11;       % 非线性环节期望损失的相角
parameter.ratiodiv = 1.025;   % 剪切频率搜寻尺度   采用* 
parameter.phidiv = 0.5;      % 相位搜寻尺度 0.5 度  采用+

parameter.later_phi = 3;     % 迟后环节相位搜寻约束系数
parameter.laterKmin = 0.5;    % 迟后搜索增益最小
parameter.laterfremin = 0.001;    % 迟后搜索中心频率最小

parameter.trapTmin = 3;
parameter.trapfremin = parameter.bandwidth * 0.2;       % 最小中心频率
parameter.trapfremax = parameter.bandwidth + 3 * pi;    % 非线性环节优化 最大频率
parameter.trapohimin = 40;  % 非线性最小相位

%% 前馈环节
parameter.para_aux1 = 1000;
parameter.para_aux2 = 1000;
parameter.forwardKmax = 0.7;


parameter.tradphi = 135;
parameter.phi_advance = 30;
parameter.phi_advance_margin = 0;
parameter.Tratio = 10;
