filename = 'ET20504162.csv';
nlength = 19;
type = 0;
% frequence = linspace(1, nlength, nlength) * 2 * pi;
if isempty(strfind(filename, '.csv')) == 0 
    ET = load(filename);
    frequence = ET(1: nlength, 1) * 2 * pi;
    mag = ET(1: nlength, 2);
    phi = ET(1: nlength, 3);
elseif isempty(strfind(filename, '.tx')) == 0 
    fid=fopen(filename,'r');
    [x, count]=fread(fid, inf, 'double'); 
    fclose(fid);
    for j = 1 : count / 3 
        frequence(j) = x(3*(j-1)+1) * 2 * pi;	
        mag(j) = x(3*(j-1)+2); 
        phi(j) = x(3*(j-1)+3);	
    end
end
response = mag.*exp(1j*phi*pi/180);
fr_data = idfrd(response,frequence,0);
if type == 0
    opt = procestOptions('Focus','simulation','Display','on','SearchMethod', 'fmincon');
else
    opt = procestOptions('Focus','simulation','Display','on','SearchMethod', 'auto');
end
sysP1D_noise = procest(fr_data,'P2I',opt);
G = tf(sysP1D_noise.Kp, [sysP1D_noise.Tp1 * sysP1D_noise.Tp2, sysP1D_noise.Tp1 + sysP1D_noise.Tp2, 1, 0]);
Mag_real = 20 * log10(mag);
Phi_real = phi;
Mag_fit = zeros(length(frequence), 1);
Phi_fit = zeros(length(frequence), 1);
[mag, phi] = bode(G, frequence);
for i = 1 : length(frequence)
    Mag_fit(i) = 20 * log10(mag(1, 1, i));
    Phi_fit(i) = phi(i);
end
figurename('对象对比2');
subplot 211;
semilogx(frequence, Mag_real, 'r*-');
hold on
semilogx(frequence, Mag_fit, 'g*-');
grid on
legend('实际', '拟合');

subplot 212;
semilogx(frequence, Phi_real, 'r*-');
hold on
semilogx(frequence, Phi_fit, 'g*-');
grid on

 
K = sysP1D_noise.Kp;
taue = sysP1D_noise.Tp1;
taum = sysP1D_noise.Tp2;
if taue > taum
	taue = taum;
	taum = sysP1D_noise.Tp1; 
end

autoArrangeFigures;