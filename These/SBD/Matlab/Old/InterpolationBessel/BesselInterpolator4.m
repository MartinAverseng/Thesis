function [ alpha,zs ] = BesselInterpolator4( derivatives,c,x)

derivatives = derivatives(:);
N = length(derivatives);
zs = BesselZeros(c,0,N,0);
A = zeros(N);
for i = 1:N
    for j = 1:N
        A(i,j) = zs(j)^(i-1)*dkBess(i-1,zs(j)*x);
    end
end

alpha = A^(-1)*derivatives;



end

