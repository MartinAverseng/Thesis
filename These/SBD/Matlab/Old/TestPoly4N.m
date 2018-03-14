%% Tester l'idée du polynome de degré 4n qui prolonge la fonction 
clear all;
close all;
N = 10;
ks = 0:2*N;
rho = 0.2;
aks = (-1).^(2*N + ks)./(rho).^(2*N + ks + 1).*cumprod(2*N+ks)/(2*N);
aks = aks(:);
x = linspace(0,rho,1000);
approx = 0;
for i = 1:length(x)
    approx(i) = sum(x(i)^(2*N)*(aks.*(x(i)-rho).^(ks'))./factorial(ks(:)));
end
x = linspace(rho/2,5*rho,1000);
plot(x,1./x);
hold on
x = linspace(0,rho,1000);
plot(x,approx);