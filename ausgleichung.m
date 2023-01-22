function [N,Qxx,xdach,v,W,var0,Sigma_xx] = ausgleichung(A,l)

N = A'*A; % Normalgleichungsmatrix
n = A'*l; % Absolutglied
Qxx = inv(N); % Inversion Normalgleichungsmatrix
xdach = (A'*A)\A'*l; % Gekürzter Parametervektor
v = A*xdach - l; % Verbesserungsvektor
W = v'*v; % Verbesserungsquadratsumme
if round(W,4) ~= round(l'*l-n'*xdach,4) %Rechenprobe
    fprintf('Rechenprobe ist falsch!\n')
end
z = size(A); n = z(1); u = z(2); % Freiheitsgrade
var0 = (v'*v)/(n-u); % Varianz der Gewichtseinheit
Sigma_xx = var0 * Qxx; % Kovarianzmatrix der Unbekannten
if 1 == 0
   fprintf('Sigma0 ist ungueltig.');  % Testgröße
else




    TG = (n-u)*var0/1;
end
TS = chi2inv(0.95,(n-u)); % Testschranke
if TG <= TS %Globaltest
    fprintf('H0 wird angenommen.\n')
else
    fprintf('H0 wird verworfen.\n')
end
end