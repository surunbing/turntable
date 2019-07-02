clc, clear
close all

load('dir_data.mat');

K = 3.491260430644539e+02;     %% ET205����
taue = 0.001726113287549;
taum = 2.465408813244642;

G_model = tf(K, [taue * taum * 1 taue + taum * 1 1 0]);

load('tra_controller.mat');
load('dir_controller.mat');

figurename('�����Ա�');
bode(G_model * K);
hold on
grid on
bode(G_model * P_trap);
legend('����ģ�͵ķ���', '���ھ���ķ���');



figurename('�ջ��Ա�');
bode(G_model * K / (1 + G_model * K));
hold on
grid on
bode(G_model * P_trap / (1 + G_model * P_trap));
legend('����ģ�͵ķ���', '���ھ���ķ���');


figurename('����');
bode(1 / (1 + G_model * K));
hold on
grid on
bode(1 / (1 + G_model * P_trap));
legend('����ģ�͵ķ���', '���ھ���ķ���');

figurename('��������');
bode(K);
hold on
grid on
bode(P_trap);
legend('����ģ�͵ķ���', '���ھ���ķ���');

figurename('��Ծ');
step(G_model * K / (1 + G_model * K), 2);
hold on
grid on
step(G_model * P_trap / (1 + G_model * P_trap), 2);
legend('����ģ�͵ķ���', '���ھ���ķ���');

autoArrangeFigures;