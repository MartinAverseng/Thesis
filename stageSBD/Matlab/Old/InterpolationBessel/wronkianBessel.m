function [ wx,A ] = wronkianBessel(N,x)


zs = BesselZeros(0,0,N,0);
A = zeros(N);
for i = 1:N
    for j = 1:N
        A(i,j) = zs(j)^i*dkBess(i,zs(j)*x);
    end
end

wx = det(A);


end

