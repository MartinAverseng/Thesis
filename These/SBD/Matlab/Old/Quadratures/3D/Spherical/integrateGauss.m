function [I] = integrateGauss(func,N)

[x,y,z,w] = Gauss_LegendreSpherical(N);
W = repmat(w,1,2*N)/(2*N);

I = sum(sum(W.*func(x,y,z)));


end

