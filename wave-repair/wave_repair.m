clc, clear
close all

%% ��ȡ�����Ƶ�����
[P, G, para] = direct_design();

%% ���������������ͺ󻷽ڽ��ͼ���Ƶ��
%% ��鲨�Σ������ݲ��˲����������
%% ����ٺ󻷽ڿ�������Ƶ�����ڽϸߵĲ�����һ����г�񣬱ջ�ʱ���ݲ��˲�������ķ�Ƶ��ʧ���Ե�����Χ��
%% �������Ƶ�������ʧ

% �������֮�󻷽�, ����wc�Ƿ����������ַ���
wc_r = 270;

wc_max = fsolve(@(x)myfun(x, para.xi, para.omegan, para.T), 100);
if wc_max > wc_r
%     wc_r = wc_r;
else
    wc_r = wc_max;
end
[mag, phi] = bode(P * G, wc_r);
phi_n_min = -123 - phi;
phi_n_max = -120 - phi;

%%�������ͺ����Ƶ����
frequence = [wc_r, wc_r - 50, para.dt];
[mag, phi] = bode(P * G, frequence);
c_data = zeros(length(frequence), 1);
for i = 1 : length(frequence)
    c_data(i) = mag(1, 1, i) * complex(cos(phi(1, 1, i) / 180 * pi), sin(phi(1, 1, i) / 180 * pi));
end


tic
% �õ��ٺ󻷽ڼ��� 
start = [wc_r / 2, 6, 2];
lb = [1; 1.5; 0.5];
ub = [wc_r; 15; 5];
% options = optimset('Algorithm','interior-point');
options = optimset('Algorithm','sqp');
[X, fval, exitflag] = fmincon(@(x)GetAlphacost(x, c_data, frequence)...
    , start, [], [], [], [], lb, ub, @(x)nonlcon(x, wc_r - 20, phi_n_min, phi_n_max, c_data, frequence), options);
toc

[c, ceq] = nonlcon(X, wc_r - 20, phi_n_min, phi_n_max, c_data, frequence);

fre = X(1);
alpha = X(2);

tau = 1 / (sqrt(alpha) * fre);
T = alpha * tau;
G_later = tf([tau, 1], [T, 1]);

figurename('��������');
K = P * G * G_later * X(3);


%% �����ݲ�����
e = 1.5;
T = 15;
f1 = 70;

trap = tf([1, e * T, f1 * f1], [1, T, f1 * f1]);
K = K * trap;

e = 2.0;
T = 15;
f1 = 88;
trap = tf([1, e * T, f1 * f1], [1, T, f1 * f1]);
K = K * trap;

e = 1.8;
T = 15;
f1 = 113;
trap = tf([1, e * T, f1 * f1], [1, T, f1 * f1]);

%% �ٺ󻷽ڣ����ӿ�������
alpha = 2;
fre = 0.5;
tau = 1 / sqrt(alpha) / fre;
T = alpha * tau;
G_later = tf([tau, 1], [T, 1]) * alpha;


K = K * trap * G_later * G_later;

margin(K);
grid on
figurename('�ջ�����');
bode(K / (1 + K));
grid on

autoArrangeFigures;
% close all
% p = tf([tau, 1], [(alpha + 1) * tau, 2]);
% bode(p);
% grid on


