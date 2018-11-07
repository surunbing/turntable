function [data, data_con] = UpdateOptdata(G_P, bandwidth, ratio, nRp, nCon, nCRp)
%   更新优化所需的data结构体
count = round(bandwidth / 2 / pi);
[fre_Rp, fre_con] = GetFrequence(bandwidth, ratio, nRp, nCon, nCRp);
fre = [linspace(1, count, count)' * 2 * pi; fre_Rp];
[mag, phi, ~] = bode(G_P, fre);

data.fre = fre;
data.mag = zeros(length(fre), 1);
data.phi = zeros(length(fre), 1);
data.phi_rad = zeros(length(fre), 1);  
data.complex = zeros(length(fre), 1);

for i = 1:length(fre)
    data.mag(i) = 20 * log10(mag(1, 1, i));
    data.phi(i) = phi(1, 1, i);
    data.phi_rad(i) = phi(1, 1, i) / 180 * pi;
    data.complex(i) = 10 ^ (mag(1, 1, i) / 20) * complex(cos(data.phi_rad(i)), sin(data.phi_rad(i)));
end

%% 约束对应的值
data_con.fre = fre_con;
data_con.mag = zeros(length(fre_con), 1);
data_con.phi = zeros(length(fre_con), 1);
data_con.phi_rad = zeros(length(fre_con), 1);  
data_con.complex = zeros(length(fre_con), 1);
[mag, phi, ~] = bode(G_P, data_con.fre);

for i = 1:length(fre_con)
    data_con.mag(i) = 20 * log10(mag(1, 1, i));
    data_con.phi(i) = phi(1, 1, i);
    data_con.phi_rad(i) = phi(1, 1, i) / 180 * pi;
    data_con.complex(i) = 10 ^ (mag(1, 1, i) / 20) * complex(cos(data_con.phi_rad(i)), sin(data_con.phi_rad(i)));
end

end

