function[out] = fromChebSeries(fun,N,X)

out = 0*X;


for n = 1:N
    out = out + fun(n)*chebyshevT(n,X);
end


end