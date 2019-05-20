function [bmag, bphi] = Getbindex(data, maglim, philim)
%% ����db ʹ����λ
bmag = 1;
bphi = 1;
for i = 1 : length(data.fre)
    if data.mag(i) >= (1 + maglim) || data.mag(i) <= (1 - maglim)
        bmag = 0;
    end
    if data.phi(i) >= philim || data.phi(i) <= -philim
        bphi = 0;   
    end
end

end

