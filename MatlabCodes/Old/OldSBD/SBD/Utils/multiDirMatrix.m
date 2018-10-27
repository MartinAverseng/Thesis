function [ M ] = multiDirMatrix( n)

A = (2*(1:n)).^2;
M = A;
for i = 1:(n-1)
    A = [0,A(1:end-1)];
    M = [M;M(end,:).*A];
end

M = [ones(1,n+1);[zeros(n,1) M]];
disp(M);


end

