w = 120 * 2 * pi;
wm = 40 * 2 * pi;

phi0 = 74 / 180 * pi;

Cost = zeros(10, 1);

for n = 2 : 10
a  = (1 - sin(phi0 / n)) / (1 + sin(phi0 / n));

cost = sqrt(w * w + a * wm * wm) / sqrt(a * a * w * w + a * wm * wm);
cost = 20 * log10(cost);
cost = cost + n;

Cost(n) = cost;


end

Cost