% Décomposition de Y0 sur les J0
clear all;
close all;
Dim = 2;
k = Dim/2 -1;
N = 200;
mu = 300;
c1 =  mu*besselj(1,mu)/besselj(0,mu);
c2 =  mu*bessely(1,mu)/bessely(0,mu);
rho1 = BesselZeros(c1,k,N,0);
rho2 = BesselZeros(c2,k,N,0);
x = linspace(0,1,100000);
Cn1 = zeros(N,1);
Cn2 = zeros(N,1);
for i = 1:N
    Cn1(i) = 1/ScalBess(rho1(i),rho1(i),k,0.00001,0.999999);
    Cn2(i) = 1/ScalBess(rho2(i),rho2(i),k,0.00001,0.999999);
end

scalYJ1 = (2/pi + mu*besselj(0,rho1)*bessely(1,mu)-rho1.*besselj(1,rho1)*bessely(0,mu))./(mu^2-rho1.^2);
scalYJ2 = (2/pi + mu*besselj(0,rho2)*bessely(1,mu)-rho2.*besselj(1,rho2)*bessely(0,mu))./(mu^2-rho2.^2);

scalYJ1(isinf(scalYJ1)) = 0; % Dans le cas où les fréquences sont égales, le produit scalaire vaut 0.
scalYJ2(isinf(scalYJ2)) = 0;

decomp1 = Cn1.*scalYJ1;
decomp2 = Cn2.*scalYJ2;

coeff1 = (Cn1.*scalYJ1)*ones(size(x));
coeff2 = (Cn2.*scalYJ2)*ones(size(x));
approx1 = sum(besselj(0,rho1*x).*coeff1,1);
approx2 = sum(besselj(0,rho2*x).*coeff2,1);
figure

plot(x,log(abs(approx2-bessely(0,mu*x))));
hold on
plot(x,log(abs(approx1-bessely(0,mu*x))));
legend({'Condition ajustée','Condition non ajustée'})
xlabel('x')
ylabel('log((Y_0 - S_N)^2)')


figure; 
plot(x,bessely(0,mu*x));
hold on;
plot(x,approx1);
plot(x,approx2);

figure;
scatter(rho1,decomp1)
hold on
scatter(rho2,decomp2)
