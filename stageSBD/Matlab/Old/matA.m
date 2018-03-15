function [ A,ddet ] = matA(N,ep)

A = zeros(N);
for i = 1:N
    for j = 1:N
        if 2*j-i -1 >= 0
            A(i,j) = factorial(2*j)/factorial(2*j-i -1)*ep^(2*j-i -1);
        end
    end
end

ddet = prod(1:2:(2*N-3))*det(A);


end

