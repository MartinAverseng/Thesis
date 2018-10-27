function [ b ] = Y0Coeffs(a,b,rho,k)

C = 1./(sqrt(pi)*rho.*abs(besselj(1,rho)));
F = @(r,rho)(2*pi*r*rho*k./(rho.^2 - k^2).*...
    ( k*bessely(0,k*r)*besselj(1,rho*r) - rho.*besselj(0,rho*r)*bessely(1,k*r) ));

b = C.*(F(b,rho) - F(a,rho));

end

