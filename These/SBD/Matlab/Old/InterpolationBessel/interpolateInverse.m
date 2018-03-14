%% Interpolate 1/x with bessel functions
addpath(genpath(pwd));
clear all
close all
N = 3; % Ordre du raccord en nombre de dérivées
rho = 0.05;
c = Inf;
derivatives = ((-1).^(0:N)).*factorial(0:N)./(rho.^(1:N+1));
derivatives(1) = 1/rho-1;
x = ChebyshevRoots(N, 'Tn',[rho 1-rho] );
x = sort(x);
pointwise = 1./x-1;

% [alpha1,zs1] = besselInterpolator(derivatives,c,rho);
[alpha3,zs3] = BesselInterpolator3(pointwise,c,x);
[alpha4,zs4] = BesselInterpolator4(derivatives,c,rho);

t = linspace(0.01,1,1000);
t3 = linspace(0,1,1000);

t4 = linspace(0,1.1*rho,1000);
% approx1 = zeros(size(t));
approx3 = zeros(size(t3));
approx4 = zeros(size(t4));
% for k = 1:length(alpha1)
%     approx1 = approx1 + alpha1(k)*besselj(0,zs1(k)*t);
%     
% end
for k = 1:length(alpha3)
    approx3 = approx3 + alpha3(k)*besselj(0,zs3(k)*t3);
    
    
end
for k = 1:length(alpha4)
    approx4 = approx4 + alpha4(k)*besselj(0,zs4(k)*t4);
    
    
end


% Afficher la fonciton et l'approximation
% plot(t,approx1);

figure;
plot(t,1./t-1);
hold on;
plot(t3,approx3);
plot(t4,approx4);

% Afficher le spectre 
figure
% scatter(zs1,log(abs(alpha1)));
% hold on;
scatter(zs3,log(abs(alpha3)));
hold on
scatter(zs4,log(abs(alpha4)));
title('Spectre (dB)')

figure
plot(t3,log(abs(1./t3 - 1 - approx3)))
hold on
plot(t4,log(abs(1./t4 - 1 - approx4)))

