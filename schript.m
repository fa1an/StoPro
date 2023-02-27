 % Stochastische Prozesse Beleg 4 - ﻿Zeitreihenanalyse von GNSS-Daten mit power-law noise
close all; clear; clc; format long

read = table2array(readtable('aufgabe2.dat'));
u = read(:,3);
t = read(:,2);
I = length(u);
if mod(I,2) == 1
    u = u(1:end-1);
    I = I-1;
end
clearvars dat

% Aufgabe 1: Sicherstellen, dass u(t) eine gerade Anzahl an Werten enthält

if mod(length(u),2)==0
else
u = u(1:end-1); % Kürzen der Zeitreihe um einen Wert
end
if mod(length(t),2)==0
else
t = t(1:end-1); % Kürzen der Zeitreihe um einen Wert
end

figure(1)
hold on
plot(t,u,'b.-')
title("GNSS Up-Reihe")
xlabel("Zeit [t]")
ylabel("u[t]")
hold off

% Aufgabe 2: Bestimmen der Kenngrößen

I = length(u); % Länge der Zeitreihe
dt = t(2)-t(1);
T = I*dt; % Abtastweite im Frequenzbereich
fg = pi/dt; % Grenzfrequenz oder Nyquistfrequenz

disp('Länge der Zeitreihe:')
disp(I);
disp('Abtastweite:')
disp(T)
disp('Grenzfrequenz:')
disp(fg)
disp('Frequenzachse:')

% Aufgabe 3: Erste Parameterschätzung

L = u;

% Naeherungsvektor der Unbekannten:
u0_ = 0;
u1_ = 1;
a1_ = 1;
a2_ = 1;
b1_ = 1;
b2_ = 1;

X0 = [u0_, u1_, a1_, b1_, a2_, b2_];

A = [];
    A(:,1) = ones(length(u),1);
    A(:,2) = t;
    A(:,3) = cos(2*pi.*t/365.25);
    A(:,4) = sin(2*pi.*t/365.25);
    A(:,5) = cos(4*pi.*t/365.25);
    A(:,6) = sin(4*pi.*t/365.25);


Li = ones(1000,1)*9999;
Ldach = ones(1000,1)*-9999;
while round(Li,6) ~= round(Ldach,6)

L0 = [];
L0(:,1) = u0_ + u1_.*t + a1_*cos(2*pi.*t/365.25) + a2_*sin(2*pi.*t/365.25) + b1_*cos(4*pi.*t/365.25) + b2_*sin(4*pi.*t/365.25);

% gekuerzter Beobachtungsvektor:
l = L - L0;

% P = Einheitsmatrix

% Ausgleichung
% ------------

[N,Qxx,xdach,v,W,var0,Sigma_xx] = ausgleichung(A,l);
s0ausg = sqrt(var0);
% Ausgeglichene Beobachtungen
Ldach = L + v;

% Hauptrechenprobe
Xdach = [ X0 ]' + xdach; % Ausgeglichene Parameter


u0_ = Xdach(1);
u1_ = Xdach(2);
a1_ = Xdach(3);
a2_ = Xdach(4);
b1_ = Xdach(5);
b2_ = Xdach(6);

Li = [];
Li(:,1) = u0_ + u1_.*t + a1_*cos(2*pi.*t/365.25) + a2_*sin(2*pi.*t/365.25) + b1_*cos(4*pi.*t/365.25) + b2_*sin(4*pi.*t/365.25);

if round(Li,6) == round(Ldach,6)
    fprintf('Hauptrechenprobe war erfolgreich.\n')
else
    fprintf('Hauptrechenprobe ist fehlgeschlagen.\n')
end

end

figure
hold on
plot(t,u,'b.-')
plot(t,Li,'r-',LineWidth=3)
title("erste Parameterschätzung")
xlabel("Zeit [t]")
ylabel("u[t]")
legend("GNSS Reihe","Parameterschätzung (mm*a^-1)","Location","southwest")
hold off

% Aufgabe 4: Darstellen der Zeitreihe der Residuen

e = v; % Berechnen der Residuen (Verbesserungen aus Ausgleichung)


figure
plot(t,e,'b.-')
title('Residuen')
xlabel('Zeit, t (d)')
ylabel('Residuen, e(t) (mm)')


% Aufgabe 5: Schätzen der Energiedichte P(f)

clearvars -except e u t I dt T fg Xdach A v Li s0ausg

E = fft(e);
E = E.*conj(E)/I;

P(1) = E(1);
P(2:I/2) = 2*E(2:I/2);
P(I/2+1) = E(I/2+1);

domega = 2*pi/T;
n = 1:round(fg/domega);
f(:,1) = n .* domega; 
f_deg = f * 180/pi;



% Aufgabe 6 - Bestimmung Spektralindex

P0 = P(1);
P_t = log(P(2:end) );
f_t = log(f);

p = polyfit(f_t,P_t,1);

fprintf("Kappa: %.4f\n", p(1))
fprintf("Kappa0: %.4f\n", p(2))

p = polyfit(log10(f),log10(P(2:end)),1);

figure
loglog(f,P(2:end),'.')
xlabel('Frequenz, f (d-1)')
ylabel('Energiedichte, P(f) (mm2)')
grid on
axis square
hold on
loglog(f,10.^(p(2)+p(1).*log10(f)),'r-')
title('Energiedichte P(f) und lineare Regressionsfunktion')
legend('P(f)','Regressionsfunktion')


% Aufgabe 7 - Kovarianzmatrix des power-law noise

U = eye(I);
h1 = 1;
for i = 2:I
    h = (h1/(i-1)) * (i - (p(1)/2) - 2);
    for j = 1:(I-i+1)
        U(j,j+i-1) = h;
    end
    h1 = h;
end
c0 = 1;
CPL = c0 * U' * U;


figure
imagesc(CPL)
colorbar
title('Spektralindex')
xlabel('Zeilenindex')
ylabel('Spaltenindex')
% Aufgabe 8 - Wiederholung Parameterschätzung

C_PL = CPL;
A = [ones(I,1) (1:I)'*dt cos(2*pi*(1:I)'*dt/365.25) sin(2*pi*(1:I)'*dt/365.25) cos(2*pi*(1:I)'*dt/(365.25/2)) sin(2*pi*(1:I)'*dt/(365.25/2))];
x = A\u;
s0 = std(u-A*x);
C = eye(I)*s0^2 + C_PL;
x = A\u;
x_cov = inv(A'*inv(C)*A)*A'*inv(C)*u;
x_cov = [x_cov(1:2); x_cov(3:4)*365.25/(2*pi); x_cov(5:6)*365.25/(4*pi)];

disp('Die Werte der geschätzten Parameter und ihre Genauigkeiten sowie die Standardabweichung der Gewichtseinheit sind:')
disp(x_cov)
disp('Standardabweichung der Geschätzten Parameter')
disp(s0)

C_x = inv(A'*inv(C)*A);
C_x_ = corrcov(C_x);
C_x = C_x./(C_x(1,1)*C_x(2,2)*C_x(3,3)*C_x(4,4)*C_x(5,5)*C_x(6,6))^(1/6);


figure
imagesc(C_x_)
title('Korrelationsmatrix der geschätzten Parameter')
xlabel('Zeilenindex')
ylabel('Spaltenindex')
colorbar

u_fit = A*x_cov;

figure
plot(t,u,'b-')
hold on
plot(t,u_fit,'r-')
title('Zeitreihe u(t) und Trajektorie des geschätzten Modells')
xlabel('Zeit, t (d)')
ylabel('Höhe, u(t) (mm)')
legend('u(t)','Modell')

disp('Vergleich der Parameter unter 3. und 8.:')
disp('Parameter 3. 8.')
disp([x x_cov])
disp('Standardabweichung aus 3.')
disp(diag(s0ausg))


%Die Parameter des Modells unter 8. sind kleiner als die Parameter unter 
% 3., da sie die Kovarianzmatrix des power-law noise berücksichtigt. 
% Die Standardabweichung der Gewichtseinheit ist unter 8. kleiner als 
% unter 3., da die Kovarianzmatrix des power-law noise die Residuen 
% besser erklärt.
