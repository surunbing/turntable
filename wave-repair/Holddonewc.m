function [later, fval, exitflag] = Holddonewc(P, G, data, wc_r, phi_dist)
global parameter
%% ��ȡ�����Ƶ�����

%% ���������������ͺ󻷽ڽ��ͼ���Ƶ��
%% ��鲨�Σ������ݲ��˲����������
%% ����ٺ󻷽ڿ�������Ƶ�����ڽϸߵĲ�����һ����г�񣬱ջ�ʱ���ݲ��˲�������ķ�Ƶ��ʧ���Ե�����Χ��
%% �������Ƶ�������ʧ
data_out = translate_data(data, P);
[~, phi] = bode_get(data_out, wc_r);
phi_n_min = -phi_dist - phi;
phi_n_max = -phi_dist + parameter.later_phi - phi;

%%�������ͺ����Ƶ����
frequence = [wc_r, wc_r * 0.7, data];%para.dt];
[mag, phi] = bode_get(data_out, frequence);
c_data = zeros(length(frequence), 1);
for i = 1 : length(frequence)
    c_data(i) = mag(1, 1, i) * complex(cos(phi(1, 1, i) / 180 * pi), sin(phi(1, 1, i) / 180 * pi));
end

tic
% �õ��ٺ󻷽ڼ��� 
start = [wc_r / 2, 6, 2];
lb = [parameter.laterfremin; 0.1; 0.1];%parameter.laterKmin
ub = [wc_r; 15; 5];
% options = optimset('Algorithm','interior-point');
options = optimset('Algorithm','sqp');
[X, fval, exitflag] = fmincon(@(x)GetAlphacost(x, c_data, frequence)...
    , start, [], [], [], [], lb, ub, @(x)nonlcon(x, wc_r, phi_n_min, phi_n_max, c_data, frequence), options);
toc

% [c, ceq] = nonlcon(X, wc_r - 20, phi_n_min, phi_n_max, c_data, frequence);

fre = X(1);
alpha = X(2);

tau = 1 / (sqrt(alpha) * fre);
T = alpha * tau;
later.G = tf([tau, 1], [T, 1]) * X(3);
later.K = X(3);
later.fre = X(1);
later.alpha = X(2);


