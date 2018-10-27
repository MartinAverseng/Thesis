Main

a = 1e-2;
rho = besselJroots(0,100);
% Méthode pour obtenir le minimum de coeff

A = d_gramMatrix(a,rho);
D = sqrt(diag(1./(diag(A))));
b = D*laplaceD2coeffs(a,rho);
M = D*A*D;
beta = M\b;
alpha1 = D*beta;
 
r = linspace(0,1,10000);
val = 0*r;

for p = 1:length(rho)
    val = val - alpha1(p)*rho(p)*Cp(rho(p))*besselj(1,rho(p)*r);
end

val = val + sum(alpha1.*rho.*Cp(rho).*besselj(1,rho*1));

loglog(r,abs(val-(1./r - 1)));

% Méthode de la dérivation terme à terme. 

A = gramMatrix(a,rho);
b = laplaceSP(a,rho);
alpha2 = A\b;

r = linspace(0,1,10000);
val = 0*r;

for p = 1:length(rho)
    val = val - alpha2(p)*rho(p)*Cp(rho(p))*besselj(1,rho(p)*r);
end

hold on
loglog(r,abs(val-1./r));
