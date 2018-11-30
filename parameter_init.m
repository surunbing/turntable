bandwidth = 18; % ��������, hz ������������
ratio = 3;

global parameter

%% ģ�Ͳ���

% K = 434;
% taum = 0.67;
% taue = 0.0035;

% K = 1.56 * 180 / pi;
% taue = 0.0039035;
% taum = 0.984871194396488;

parameter.K = 496.7296;
parameter.taum = 2.0624;
parameter.taue = 0.0019;

parameter.bandwidth = bandwidth * 2 * pi;
parameter.ratio = ratio;        % ѡ��ԽƵ�ʱ���
parameter.wc_max = 650;     % ֱ����Ƶ�ʱ�����������ƣ������޸�
parameter.T = 0.0014;       % ֱ�����T
parameter.Tmax = parameter.T / 1.6; % ֱ����Ƶ����T 
parameter.Tdiv = 1.05; %ֱ���漰����ϵ����������
parameter.kgr = 5;     % ֱ���漰��λԣ������
parameter.Mre = 6;     % ֱ�����г�������
parameter.pmr = 35;    % ֱ�������λԣ������
parameter.dampmax = 1; % ֱ�������������� 
parameter.pmmax = 45;

%% �������͵���ز���
parameter.phi_creg = 8;   % �ջ�������λĿ��
parameter.mag_creg = 0.8; % �ջ����η�ֵĿ��, db
parameter.maglim = 0.45;  % �ջ����շ�ֵ db
parameter.philim = 5;     % �ջ�������λ
parameter.num_max = 3;    % �����������Ի�������
parameter.rdiv = 1.05;    % ָ���Ż���ſ���ϵ��
parameter.start_ratio = parameter.ratio * 0.8;
parameter.phimarginmin = 47.5;  % ��λԣ����Сֵ

parameter.phi_margin = 122;  %�ջ�����Ϊ��ʹ���ղ�����, ��Ҫ�����Ż�����λԣ��    �����ֵ����Ѱ�� ������ʵ�飬Ч������
parameter.phi_reg = 8;       % �����Ի���������ʧ�����
parameter.ratiodiv = 1.025;   % ����Ƶ����Ѱ�߶�   ����* 
parameter.phidiv = 0.5;      % ��λ��Ѱ�߶� 0.5 ��  ����+

parameter.later_phi = 3;     % �ٺ󻷽���λ��ѰԼ��ϵ��
parameter.laterKmin = 0.5;    % �ٺ�����������С
parameter.laterfremin = 0.001;    % �ٺ���������Ƶ����С

parameter.trapTmin = 3;
parameter.trapfremin = parameter.bandwidth * 0.2;       % ��С����
parameter.trapfremax = parameter.bandwidth + 3 * pi;    % �����Ի����Ż� ���Ƶ��
parameter.trapohimin = 40;  % ��������С��λ

%% ǰ������
parameter.para_aux1 = 500;
parameter.para_aux2 = 500;
parameter.forwardKmax = 0.7;
