bandwidth = 10; % ��������, hz ������������
ratio = 3.0;

global parameter

%% ģ�Ͳ���
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
% outputfile = 'C:\Users\Momenta\Documents\��ҵ���\����\CSDA_FANGXUN\turntable/controller5.txt';

parameter.K = K;
parameter.taum = taum;
parameter.taue = taue;

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
parameter.maglim = 0.7;  % �ջ����շ�ֵ db
parameter.philim = 8;     % �ջ�������λ
parameter.num_max = 3;    % �����������Ի�������
parameter.rdiv = 1.05;    % ָ���Ż���ſ���ϵ��
parameter.start_ratio = parameter.ratio * 0.8;
parameter.phimarginmin = 44.5;  % ��λԣ����Сֵ

parameter.phi_margin = 122;  %�ջ�����Ϊ��ʹ���ղ�����, ��Ҫ�����Ż�����λԣ��    �����ֵ����Ѱ�� ������ʵ�飬Ч������
parameter.phi_reg = 11;       % �����Ի���������ʧ�����
parameter.ratiodiv = 1.025;   % ����Ƶ����Ѱ�߶�   ����* 
parameter.phidiv = 0.5;      % ��λ��Ѱ�߶� 0.5 ��  ����+

parameter.later_phi = 3;     % �ٺ󻷽���λ��ѰԼ��ϵ��
parameter.laterKmin = 0.5;    % �ٺ�����������С
parameter.laterfremin = 0.001;    % �ٺ���������Ƶ����С

parameter.trapTmin = 3;
parameter.trapfremin = parameter.bandwidth * 0.2;       % ��С����Ƶ��
parameter.trapfremax = parameter.bandwidth + 3 * pi;    % �����Ի����Ż� ���Ƶ��
parameter.trapohimin = 40;  % ��������С��λ

%% ǰ������
parameter.para_aux1 = 1000;
parameter.para_aux2 = 1000;
parameter.forwardKmax = 0.7;


parameter.tradphi = 135;
parameter.phi_advance = 30;
parameter.phi_advance_margin = 0;
parameter.Tratio = 10;
