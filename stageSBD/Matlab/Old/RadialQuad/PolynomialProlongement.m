%k;% Approximation de la fonction inverse selon la stratégie du prolongement polynomial :
clear all;
close all;

Dim = 2;
N = 6;
a = 0.05;
b = 1-a;

%f = 1/x;
ak = zeros(2*N+1,1);
for k = 0:(2*N)
    ak(k+1) = pochhammer(-k-2*N, k)*a^(-1-k-2*N);
end


bk = zeros(2*N+1,1);
for k = 0:(2*N)
    bsum = 0;
    for l = 0:k
        bsum = bsum+nchoosek(2*N+l-1,l)*(b/(b-1))^l;
    end
    bk(k+1) = (-1)^k*factorial(k)/b^(k+1)/(1-b)^(2*N)*bsum;
end

% Check graphically validity of extension

func = @(x)(1./x);
ext = @(x)(polyExt(func,a,b,ak,bk,x));
t = linspace(0,1,1000);
plot(t,ext(t));


% Compute the Bessel-Fourier series of this function
P = 250;
k = Dim/2-1;
c = Inf; % Dirichlet condition
bessZs = BesselZeros(c,k,P,0);

alpha = zeros(P,1);
for p = 1:P
    [x,w] = Gauss_Legendre1D(2000,0,1);
    zp = bessZs(p);
    e_p = besselj(k,zp*x)*sqrt(2)./(x.^k)/abs(besselj(k+1,zp));
    assert(abs(sum(w.*x.^(Dim-1).*e_p.^2)-1)<1e-7);
    alpha(p) = sum(w.*x.^(Dim-1).*ext(x).*e_p);
end

% Plot the spectrum :

scatter(bessZs,alpha);

% Compare function and approximation

approx = @(x)(besselApprox(bessZs,alpha,k,x));
t1 = linspace(a/2,1-(1-b)/2,1000);
t2 = linspace(0,1,1000);
figure
plot(t1,1./t1);
hold on
plot(t2,ext(t2));
hold on
plot(t2,approx(t2));

% Compute error
terr = linspace(a,b,1000);
eps = max(abs(1./terr-approx(terr)));
