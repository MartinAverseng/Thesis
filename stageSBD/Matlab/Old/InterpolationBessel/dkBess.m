function [res] = dkBess(k,rho)

res = 0;
for j = 0:k
    res = res + (-1)^j*nchoosek(k,j)*besselj(2*j-k,rho);
end
res = res/2^k;


end

