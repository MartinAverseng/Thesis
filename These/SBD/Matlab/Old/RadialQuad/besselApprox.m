function [ res ] = besselApprox( z,alpha,k,x,b)

if b
res = zeros(size(x));
for p = 1:length(z);
    res = res + alpha(p)*besselj(k,z(p)*x)./(x.^(k))*sqrt(2)/abs(besselj(k+1,z(p)));
end

else
res = zeros(size(x));
for p = 1:length(z);
    res = res + alpha(p)*besselj(k,z(p)*x);
end

end

