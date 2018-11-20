function [Rp] = GetRp(x, T)
%   »ñµÃĞ³Õñ·å
xi = x(1);
omegan = x(2);
Rp = 0;
if xi >= 1 / sqrt(2)
    Rp = 1;
else
    y = 1 - 2 * xi * xi;
    ar = 3 * T * T;
    br = 2 - 4 * T * T * omegan * omegan * (1 - 2 * xi * xi);
    cr = omegan * omegan * (T * T * omegan * omegan + 4 * xi * xi - 2);
    if br > 0 && cr >= 0
        Rp = 1;
    elseif br > 0 && cr < 0
        if xi > 0 && xi < 1 / 2 && T * T * omegan * omegan < (1 / (2 * y))
            Rp = -1;
        elseif xi > 1 / 2 && xi < 1 / sqrt(2) && T * T * omegan * omegan < (2 * y)
            Rp = -1;
        end
    elseif br < 0 && cr > 0
        para1 = T * T * omegan * omegan;
        para2 = (sqrt(3 * (4 * T ^ 4 * omegan ^ 4 - 1)) - 1) / (4 * para1);
        para3 = 2 * y;
        if para1 > para3 && xi > 0 && xi < 1 / 2 && y > para2
            Rp = -1;
        elseif para1 > 1 / para3 && para1 >= 1 / 2 && xi >= 1 / 2 && xi < 1 / sqrt(2) && y > para2
            Rp = -1;
        end
    elseif br < 0 && cr <= 0
         para1 = T * T * omegan * omegan;
         if para1 < 2 * y && para1 > 1 / (2 * y)
             Rp = -1;
         end
    end
    if Rp == -1
        deltar = (4 * T * T * omegan * omegan * y + 1) ^ 2 - 3 * (4 * T ^ 4 * omegan ^ 4 - 1);
        x1 = (-br + sqrt(deltar)) / (2 * ar);
        f = (omegan * omegan - x1) ^ 2 + 4 * xi * xi * omegan * omegan * x1;
        delta = (2 * T * T * omegan * omegan * (1 - 2 * xi * xi)) / 3 + 2 / 3;
        Rp = omegan * omegan / sqrt(f * (sqrt(deltar) / 6 + delta));
    else
        Rp = 1;
        
    end
end
end

