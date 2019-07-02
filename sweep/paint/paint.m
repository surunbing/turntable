clc, clear
close all

load('dir_data.mat');

K = 3.491260430644539e+02;     %% ET205半载
taue = 0.001726113287549;
taum = 2.465408813244642;

G_model = tf(K, [taue * taum * 1 taue + taum * 1 1 0]);

load('tra_controller.mat');
load('dir_controller.mat');

figurename('开环对比');
bode(G_model * K);
hold on
grid on
bode(G_model * P_trap);
legend('基于模型的方法', '基于经验的方法');



figurename('闭环对比');
bode(G_model * K / (1 + G_model * K));
hold on
grid on
bode(G_model * P_trap / (1 + G_model * P_trap));
legend('基于模型的方法', '基于经验的方法');


figurename('噪声');
bode(1 / (1 + G_model * K));
hold on
grid on
bode(1 / (1 + G_model * P_trap));
legend('基于模型的方法', '基于经验的方法');

figurename('噪声开环');
bode(K);
hold on
grid on
bode(P_trap);
legend('基于模型的方法', '基于经验的方法');

figurename('阶跃');
step(G_model * K / (1 + G_model * K), 2);
hold on
grid on
step(G_model * P_trap / (1 + G_model * P_trap), 2);
legend('基于模型的方法', '基于经验的方法');

autoArrangeFigures;