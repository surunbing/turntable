% G ∂‘œÛ
function [c,ceq] = nonlcon1(x, data, series, P, bandwidth, fre)
G = GetTf(x, series);
[mag, phi, w] = bode(G * P);
[Gm, Pm, Wcg, Wmp] = margin(G * P);
wmin = abs(w - Wmp);
num = find(wmin == min(wmin));
phi = phi(1:num);
phi_min = min(phi);
c(1) = Wmp - 320;
c(2) = 35 - Pm;
c(3) = 1 - Gm;
c(4) = Pm - 55;
c(5) = 150 - Wmp;
c(6) = abs(phi_min) - 170;
ceq = [ ];