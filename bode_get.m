function [mag, phi] = bode_get(data, fre)
% ������ɢ���ݣ���ȡָ����ķ�Ƶ����Ƶ����
mag = interp1(data.fre, data.mag, fre);
phi = interp1(data.fre, data.phi, fre);
end

