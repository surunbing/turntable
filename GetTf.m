function [G] = GetTf(x, series)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
G = x(1);

for i = 1:series.real_pole
    G = G * tf(1, [1 / x(i + 1), 1]);
end

num = series.real_pole + 1;
for i = 1:series.real_zero
    P = tf([1 / x(i + num), 1], 1); 
    G = G * P; 
end

num = num + series.real_zero;
for i = 1:series.trap
    e = x(i * 4 + num - 3);
    T = x(i * 4 + num - 2);
    f = x(i * 4 + num - 1);
    f1 = x( i * 4 + num); 
    P = tf([1, e * T, f^2] / f / f, [1, T, f1^2] / f1 / f1);
    G = G * P; 
end

num = num + series.trap * 4;
for i = 1:series.lead  
%     P = tf([1 / x(i * 2 + num - 1), 1], [1 / x(i * 2 + num), 1]);
    alpha = x(i * 2 + num - 1);
    fre = x(i * 2 + num);
    tau = 1 / (sqrt(alpha) * fre);
    T = alpha * tau;
    P = tf([tau, 1], [T, 1]);
    G = G * P; 
end


end