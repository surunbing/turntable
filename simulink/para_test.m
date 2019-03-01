wq = 150 * 2 * pi;

% Q1 = tf(0.2756 * wq ^ 4, [1, 0.9528 * wq, 1.4539 * wq * wq, 0.7426 * wq ^ 3, 0.2756 * wq ^ 4]);
% Q2 = tf(wq ^ 4, [1, 4 * wq, 6 * wq * wq, 4 * wq ^ 3, wq ^ 4]);
% bode(Q1);
% grid on
% hold on
% bode(Q2);

% J = 1.076;
% L = 0.733240176511722;
% r = 80.116666666666666667;
% Kt = 24.812626200641272;
% Kv = 0.12;
% Kps = 1.077250780790937e+02;

% Type = 1; %% sine
% t_off = 10;
% mag = 5;
% fre = 1;

K = 1.56 * 180 / pi;
taue = 0.0039035;
taum = 0.984871194396488;
b2  = 1 / taue;
b1 = 1 / taum;
b0 = K / (taue * taum);

%% «–±»—©∑Ú
beta4 = 0.2756 * wq ^ 4;
beta1 = 0.9528 * wq - b1 - b2;
beta2 = 1.4539 * wq * wq - b2 * beta1 - b1 * beta1 - b1 * b2;
beta3 = 0.7426 * wq ^ 3 - b2 * beta2 - b1 * beta2 - b1 * b2 * beta1;
Q1 = tf(beta4, [1, 0.9528 * wq, 1.4539 * wq^2, 0.7426 * wq^3, 0.2756 * wq^4]);
bode(Q1);
grid on
hold on
gain = b0;

%% º´µ„≈‰÷√
% beta4 = wq ^ 4;
% beta1 = 4 * wq - b1 - b2;
% beta2 = 6 * wq * wq - b2 * beta1 - b1 * beta1 - b1 * b2;
% beta3 = 4 * wq ^ 3 - b2 * beta2 - b1 * beta2 - b1 * b2 * beta1;
% Q2 = tf(beta4, [1, 4 * wq, 6 * wq^2, 4 * wq^3, 1 * wq^4]);
% bode(Q2);
% grid on

