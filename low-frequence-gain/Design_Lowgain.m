% load('control.mat');
%% �����Ƶ���棬ͬʱ���ƻ���Ƶ��λ
%% ��ʹ�û�����, �Ƕ��ƶ�

%% ����ṹ��
LowGain.count = 0;
LowGain.fre = 0;
LowGain.alpha = 0;
LowGain.K = 1;
LowGain.Pmloss = 0;

%% ѡ��Ƶ�㣬�����Ƶ������������λ˥�������ԣ�ȵ�Ӱ��
%% Ƶ�㲻��̫�ͣ�����Ӱ�������ľ�̬����
LOW_FRE_MIN = 0.05 * 2 * pi;
LOW_FRE_MAX = 1 * 2 * pi;
LOW_PMLOSS_MAX = 1;

%% �����Ƶ����Ӱ����λԣ�ȹ���
%% ����Ƶ��, ��С���Ƶ������

LowGain = GetNewgain(LowGain, LOW_FRE_MIN, LOW_FRE_MAX, 3, 250);
LowGain = GetNewgain(LowGain, LOW_FRE_MIN, LOW_FRE_MAX, 3, 250);
LowGain = GetNewgain(LowGain, LOW_FRE_MIN, LOW_FRE_MAX, 3, 250);


Glow = GetlowgainG(LowGain);
figurename('�ٺ󻷽�');
bode(Glow);
grid on

gain = 20 * log10(27)
