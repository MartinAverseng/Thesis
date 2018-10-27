function [ M ] = symmetrize(M)

M = M + triu(M,1)' - tril(M,-1);


end

