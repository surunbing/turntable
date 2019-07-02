clc, clear
close all

K = 3.491260430644539e+02;     %% ET205Аыди
taue = 0.001726113287549;
taum = 2.465408813244642;

G_model = tf(K, [taue * taum * 1 taue + taum * 1 1 0]);

Num = cell2mat(G_model.Numerator);
Den = cell2mat(G_model.Denominator);

load('tra_controller.mat');
load('dir_controller.mat');

KNum = cell2mat(K.Numerator);
Kden = cell2mat(K.Denominator);

PNum = cell2mat(P_trap.Numerator);
Pden = cell2mat(P_trap.Denominator);


