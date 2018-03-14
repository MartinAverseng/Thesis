function [ alpha,zs ] = besselInterpolator2(ft,fpt,c,t)
t = t(:);
ft = ft(:);
fpt = fpt(:);
zs = BesselZeros(c,0,2*length(t),0);
N = length(t);
A = zeros(2*N);
for i = 1:2*length(ft)
    for j = 1:2*length(ft)
        if i<=N
            A(i,j) = besselj(0,zs(j)*t(i));
        else
            A(i,j) = -zs(j)*besselj(1,zs(j)*t(i-N));
        end
    end
end
alpha = A\[ft;fpt];

end

