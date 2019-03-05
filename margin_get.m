function [Gm, Pm, Wgm, Wpm] = margin_get(data)
%���ݿ������ݻ�ü���Ƶ�ʣ����ԣ�ȵȵ�
Wgm = interp1(data.phi, data.fre, -180);
Gm = interp1(data.fre, data.mag, Wgm);
Wpm = interp1(data.mag, data.fre, 0);
Pm = interp1(data.fre, data.phi, Wpm);
end

