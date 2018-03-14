function [ ak ] = computeAk(n,a)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

ak = zeros(2*n+1,1);
for k = 0:2*n
    for j = 0:k-1
        ak(k+1) = ak(k+1) + nchoosek(j+2*n-1,j)/(k-j);
    end
    ak(k+1) = (-1)^k*factorial(k)/a^(2*n+k)*(-ak(k+1) + nchoosek(k+2*n-1,k)*log(a));
end


end

