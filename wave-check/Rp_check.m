function [fre, Rp, flag] = Rp_check(data, bandwidth, ratio, Rp_max, option)
%   ·µ»Ø±Õ»·Ð³Õñ·å, db
if strcmp(option.type, 'transfer') == 1
    [mag, phi, frequence] = bode(data);
    ncount = length(frequence);
    Mag = zeros(ncount, 1);
    Phi = zeros(ncount, 1);
    for i = 1 : ncount
        Mag(i) = mag(1, 1, i);
        Phi(i) = phi(1, 1, i);
    end
    Rp = max(Mag);
    num = find(Mag == Rp);
    fre = frequence(num(1));
    Rp = 20 * log10(Rp);
elseif strcmp(option.type, 'discrete') == 1
    Rp = max(data.mag);
    num = find(data.mag == Rp);
    fre = data.fre(num(1)); 
    Rp = 20 * log10(Rp);
end
if Rp > Rp_max
    flag = 0;
elseif Rp <= Rp_max
    flag = 1;
end
end

