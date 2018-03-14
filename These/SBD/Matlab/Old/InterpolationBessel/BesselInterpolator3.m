function [ alpha,zs ] = BesselInterpolator3( ft,c,t )
t = t(:);
ft = ft(:);

zs = BesselZeros(c,0,length(t),0);
A= besselj(0,t*zs');

alpha = A^(-1)*ft;

end

