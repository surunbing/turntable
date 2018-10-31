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
for i = 1:series.complex_pole
    pole1 = complex(-x(i * 2 + num - 1), x(i * 2 + num)); 
    pole2 = complex(-x(i * 2 + num - 1), - x(i * 2 + num)); 
    circle_rad = x(i * 2 + num - 1) ^ 2 + x(i * 2 + num) ^ 2;
    P = tf(1, conv([1, -pole1], [1, -pole2]) / circle_rad);
    G = G * P; 
end

num = num + series.complex_pole * 2;
for i = 1:series.complex_zero
    pole1 = complex(-x(i * 2 + num - 1), x(i * 2 + num)); 
    pole2 = complex(-x(i * 2 + num - 1), - x(i * 2 + num)); 
    circle_rad = x(i * 2 + num - 1) ^ 2 + x(i * 2 + num) ^ 2;
    P = tf(conv([1, -pole1], [1, -pole2]) / circle_rad, 1);
    G = G * P; 
end

num = num + series.complex_zero * 2;
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

