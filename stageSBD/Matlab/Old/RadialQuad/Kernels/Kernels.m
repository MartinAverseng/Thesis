%% Quadrature radiale de différents noyaux

clear all;
close all;


Dim = 2;
D = 50*(rand(1)+0.5);
Rmin = D/10; 
Rmax = D; % On prend de la marge
tol = 1e-5;
% le champ, dans R^2

% Noyau de Laplace
switch Dim
    case 2
        G1 = @(X)(-log(X));
    case 3
        G1 = @(X)(1./X);
end

% Noyau de Helmholtz (partie réelle)
k = 80/Rmax;
switch Dim
    case 2
        G2 = @(X)(bessely(0,k*X));
        Gdiff2 = @(X)(-bessely(1,k*X));
    case 3
        G2 = @(X)(cos(k*X)./X);
end

% On veut calculer gk = \sum_{j = 1..Nx}G(xk,yj)f(xk)

%% Quadrature radiale
    

% Echelle normalisée : 
rho = Rmin/Rmax;
G1_norm = @(X)(G1(Rmax*X) - G1(Rmax));


% Compute the Fourier-Bessel decomposition
c = Inf; % We use the dirichlet quadrature for Laplace
[beta1,freqs1,quad1,res1] = BesselQuadSchmidt(Dim,c,G1_norm,rho,1-rho,tol);

k = k*Rmax;
if abs(besselj(Dim/2-1,k))<1e-7
    c = Inf;
else
    if Dim/2-1==0
        derBess = @(x)(-bessely(1,x));
    else
        derBess = @(x)(0.5*(besselj(Dim/2-2,x)-besselj(Dim/2,x)));
    end
    c = -k*derBess(k)/bessely(Dim/2-1,k);
end
G2_norm = @(X)(G2(Rmax*X));

[beta2,freqs2,quad2,res2] = BesselQuadSchmidt(Dim,c,G2_norm,rho,1-rho,tol,k);
Cn = zeros(length(freqs2),1);

for i = 1:length(Cn)
    Cn(i) = 1/ScalBess(freqs2(i),freqs2(i),Dim/2-1,0.0000001,0.9999999);
end

scalYJ1 = (2/pi + k*besselj(0,freqs2)*bessely(1,k)-freqs2.*besselj(1,freqs2)*bessely(0,k))./(k^2-freqs2.^2);
hold on
scatter(freqs2,Cn.*scalYJ1);
legend({'Approx','Exact'})
title('Spectre pour Y_0')

x = linspace(0,1,100000);
coeff = (Cn.*scalYJ1)*ones(size(x));
approx = sum(besselj(0,freqs2*x).*coeff,1);
figure(3)
hold on
plot(x,approx);
legend({'Fonction exacte','Approximation SCSD','Troncature de la série'})
title('Y_0');



