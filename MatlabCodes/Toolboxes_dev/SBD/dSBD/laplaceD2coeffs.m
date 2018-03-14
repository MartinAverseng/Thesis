function [ b ] = laplaceD2coeffs( a,rho)

rho = rho(:);
P = length(rho);
if length(a) == 2
    b = a(2);
    a = a(1);
else
    b = 1;
end
C = Cp(rho);
dek= @(r)(-C.*rho.*besselj(1,rho*r));
b = 2*pi*(1/a*dek(a) - dek(1));


end

