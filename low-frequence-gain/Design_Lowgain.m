clc, clear
close all

%% �����Ƶ���棬ͬʱ���ƻ���Ƶ��λ
%% ��ʹ�û�����, �Ƕ��ƶ�

%% ����ṹ��
LowGain.count = 0;
LowGain.fre = 0;
LowGain.alpha = 0;
LowGain.K = 0;
LowGain.Sum = 0;
LowGain.Pmloss = 0;

%% ѡ��Ƶ�㣬�����Ƶ������������λ˥�������ԣ�ȵ�Ӱ��
%% Ƶ�㲻��̫�ͣ�����Ӱ�������ľ�̬����
LOW_FRE_MIN = 0.01 * 2 * pi;
LOW_FRE_MAX = 0.5 * 2 * pi;
LOW_PMLOSS_MAX = 3;

%% �����Ƶ����Ӱ����λԣ�ȹ���
%% ����Ƶ��, ��С���Ƶ������

