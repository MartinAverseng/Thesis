function [ alpha,zs ] = besselInterpolator(vec,c,rho)

vec = vec(:);
N = length(vec);
zs = BesselZeros(c,0,N,0);
A = zeros(N);
for i = 0:N-1
    for j = 0:N-1
        A(i+1,j+1) = zs(j+1)^i*dkBess(i,zs(j+1)*rho);
    end
end

alpha = A\vec;

end

