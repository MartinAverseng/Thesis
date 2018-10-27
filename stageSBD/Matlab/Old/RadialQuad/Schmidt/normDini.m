function [ H ] = normDini(c,gamma,k)



if c ==Inf
    % Returns the value of sqrt(\int_0^1 x J_k(gamma x)^2)
    % where gamma is assumed to be a root of J_k
    H = 1/sqrt(2)*abs(besselj(k+1,gamma));
else
    % Returns the value of sqrt(\int_0^1 x J_k(gamma x)^2)
    % where gamma is assumed to be a root of c J_k + x Jk'
    H = sqrt(c^2 + gamma^2 - k^2)*abs(besselj(k,gamma))/(sqrt(2)*gamma);
end

