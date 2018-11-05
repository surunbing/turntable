function [mag, phi] = GetMagPhi(x, series, fre)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
nLength = length(fre);
mag = ones(nLength, 1);
phi = zeros(nLength, 1);
for i = 1:series.real_pole
    res = 1 ./ complex(ones(nLength), (-1 / x(i + 1)) * fre);
    mag = mag .* abs(res);
    phi = phi + angle(res);
end

num = series.real_pole + 1;
for i = 1:series.real_zero
    res = complex(ones(nLength), (-1 / x(i + num)) *  fre);
    mag = mag .* abs(res);
    phi = phi + angle(res);
end

num = num + series.real_zero;
for i = 1:series.complex_pole
    pole1 = complex(-x(i * 2 + num - 1), x(i * 2 + num)); 
    pole2 = complex(-x(i * 2 + num - 1), - x(i * 2 + num)); 
    circle_rad = x(i * 2 + num - 1) ^ 2 + x(i * 2 + num) ^ 2;
    para = conv([1, -pole1], [1, -pole2]) / circle_rad;
    res = 1 ./ complex(ones(nLength, 1) - para(1) * fre .* fre, para(2) * fre);
    mag = mag .* abs(res);
    phi = phi + angle(res);   
end

num = num + series.complex_pole * 2;
for i = 1:series.complex_zero
    pole1 = complex(-x(i * 2 + num - 1), x(i * 2 + num)); 
    pole2 = complex(-x(i * 2 + num - 1), - x(i * 2 + num)); 
    circle_rad = x(i * 2 + num - 1) ^ 2 + x(i * 2 + num) ^ 2;
    para = conv([1, -pole1], [1, -pole2]) / circle_rad;
    res = complex(ones(nLength, 1) - para(1) * fre .* fre, para(2) * fre);
    mag = mag .* abs(res);
    phi = phi + angle(res);
end

num = num + series.complex_zero * 2;
for i = 1:series.lead  
    alpha = x(i * 2 + num - 1);
    frequence = x(i * 2 + num);
    tau = 1 / (sqrt(alpha) * frequence);
    T = alpha * tau;
    res = complex(ones(nLength, 1), tau * fre) ./ complex(ones(nLength, 1),  T * fre);
    mag = mag .* abs(res);
    phi = phi + angle(res);
end

mag = 20 * log10(mag * x(1));
phi = phi / pi * 180;

end

