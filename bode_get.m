function [mag, phi] = bode_get(data, fre)
% 给定离散数据，获取指定点的幅频和相频特性
mag = interp1(data.fre, data.mag, fre);
phi = interp1(data.fre, data.phi, fre);
end

