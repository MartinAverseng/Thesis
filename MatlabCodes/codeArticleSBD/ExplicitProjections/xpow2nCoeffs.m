function [ b ] = xpow2nCoeffs(a,rho,n)

F = @(x)(-2*n*(2*besselj(1, rho*x).*LommelS(2*n-1, 0, rho*x)*n-LommelS(2*n, 1, rho*x).*besselj(0, rho*x)).*rho.^(-2*n+1)*x);
b = 2*pi*Cp(rho).*(F(1) - F(a));

end

