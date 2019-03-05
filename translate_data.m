function data_out = translate_data(data, P)
% 一组数据加一个控制器变成一组数据
[mag, phi] = bode(P, data.fre);
mag = 20 * log10(reshape(mag, [length(data.fre), 1]));
phi = reshape(phi, [length(data.fre), 1]);
data_out.fre = data.fre;
data_out.mag = mag + data.mag;
data_out.phi = phi + data.phi;
end

