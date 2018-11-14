clc, clear
close all

%% ��ȡ�����Ƶ�����
[P, G, para] = direct_design();

%% ���������������ͺ󻷽ڽ��ͼ���Ƶ��
%% ��鲨�Σ������ݲ��˲����������
%% ����ٺ󻷽ڿ�������Ƶ�����ڽϸߵĲ�����һ����г�񣬱ջ�ʱ���ݲ��˲�������ķ�Ƶ��ʧ���Ե�����Χ��
%% �������Ƶ�������ʧ

% �������֮�󻷽�, ����wc�Ƿ����������ַ���
wc_r = 210;

wc_max = fsolve(@(x)myfun(x, para.xi, para.omegan, para.T), 100);
if wc_max > wc_r
%     wc_r = wc_r;
else
    wc_r = wc_max;
end
[mag, phi] = bode(P * G, wc_r);
phi_n_min = -130 - phi;
phi_n_max = -126 - phi;

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
lb = [1; 1.5; 1];
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

%�Ż�����K��ʹ�ü���Ƶ��
% count = 50;
% K = linspace(1, 5, count);
% Wc = zeros(count, 1);
% Pm = zeros(count, 1);
% Phi = zeros(count, 1);
% for i = 1 : count
%     [Gm, pm, Wgm, Wpm] = margin(K(i) * G_later * P * G);
%     Wc(i) = Wpm;
%     Pm(i) = pm;
%     if pm < 50 || pm > 55
%         Pm(i) = 0;
%     end
%     if Wpm > wc_r || Wpm < (wc_r - 50)
%         Wc(i) = 0;
%     end
%     [~, Phi(i)] = bode((K(i) * G_later * P * G) / (1 + K(i) * G_later * P * G), para.dt);
% end

% ��������K


figurename('��������');
K = P * G * G_later * X(3);
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


