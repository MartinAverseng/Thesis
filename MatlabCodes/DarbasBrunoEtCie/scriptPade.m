% Test de l'approximation de Padé de la racine avec rotation de branche en
% fonction des paramètres
clear all;
close all;

Np = 20;
xmax = 20*Np;
y = linspace(-xmax,xmax,500);
h1 = figure;
k = 20;
epsilon = k^(1/3);
reg = 1/(1+1i*epsilon/k)^2;
for theta = 0:0.1:pi

[ C0,Aj,Bj ] = rotatingPadeRacine(Np,theta);
Rp = @(x)(C0 + sum(Aj*x./(1 + Bj*x)));

figure(h1);
semilogy(y,abs(1i*k*Rp(-y.^2/(k+1i*epsilon)^2)-1i*k*sqrt(1 - y.^2./(k+1i*epsilon)^2)));
hold on

end


