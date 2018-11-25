function [cost] = GetAlphacost(x, c_data, frequence)
%	��С�����
% alpha = x(2);
% cost = alpha;
fre = abs(x(1));
alpha = abs(x(2));
K = abs(x(3));
tau = 1 / (sqrt(alpha) * fre);
cd_data = K * complex(1, tau * frequence) ./ complex(1, alpha * tau * frequence);
cost = abs(angle(cd_data(3) * c_data(3) / (1 + cd_data(3) * c_data(3))) / pi * 180);
end

