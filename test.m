clc, clear
close all
tic
G = tf(1, [0.1, 1]);
G = G * G * G;
[mag, phi, fre] = bode(G);
toc
tic
mag = zeros(1, 100);
phi = zeros(1, 100);
for i = 1 : 100
    a = 1 / complex(i, 1) * 1 / complex(i, 1) * 1 / complex(i, 1);
    mag(i) = abs(a);
    phi = angle(a);
end
toc

tic
a = conv([1, 2], [1, -2]) / 1;
toc