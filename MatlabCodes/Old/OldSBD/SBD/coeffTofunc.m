function [out] = coeffTofunc(alpha,rho,x)

P = length(rho);
out = 0*x;

for p = 1:P
    freq_p = alpha(p)*besselj(0,rho(p)*x)*Cp(rho(p));
    out = out + freq_p;
end


end

