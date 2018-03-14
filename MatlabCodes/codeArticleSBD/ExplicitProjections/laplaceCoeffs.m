function [ coeffs ] = laplaceCoeffs( a,b,rho)

coeffs = 2*pi*Cp(rho).*(besselj(0,rho*b) - besselj(0,rho*a));


end

