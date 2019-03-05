function [Gm, Pm, Wgm, Wpm] = margin_get(data)
%根据开环数据获得剪切频率，相角裕度等等
Wgm = interp1(data.phi, data.fre, -180);
Gm = interp1(data.fre, data.mag, Wgm);
Wpm = interp1(data.mag, data.fre, 0);
Pm = interp1(data.fre, data.phi, Wpm);
end

