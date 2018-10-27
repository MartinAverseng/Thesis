function [ vec ] = derVec( N,ep )
% Return vector vec(i) =  d^i f/dx (ep) where f = x -> 1/x

vec = zeros(N,1);
for i = 1:N
    vec(i) = (-1)^(i-1)*factorial(i-1)/ep^(i);
end


end

