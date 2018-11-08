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
for i = 1:series.trap
    e = x(i * 4 + num - 3);
    T = x(i * 4 + num - 2);
    f = x(i * 4 + num - 1);
    f1 = x(i * 4 + num);
    res = (complex(f * f * ones(nLength, 1) - fre .* fre, e * T * fre) / f / f) ./ (complex(f1 * f1 * ones(nLength, 1) - fre .* fre, T * fre) / f1 / f1);
    mag = mag .* abs(res);
    phi = phi + angle(res); 
end

num = num + series.trap * 4;
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