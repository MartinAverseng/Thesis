function [ b ] = invxCoeff(a,rho)

F=@(x)((1/2)*rho.^2*x.*(besselj(0, rho*x)*pi.*StruveH_moins1(rho*x)-2*besselj(1, rho*x).*LommelS(-2, 0, rho*x)));
b = 2*pi*Cp(rho).*(F(1) - F(a));

end

