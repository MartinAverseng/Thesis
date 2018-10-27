function [ out ] = plateau01(x)

out = 0*x;
out(x<=0) = 1;
idx = and(x>0,x<1);
out(idx) = 1-exp(-1./(x(idx)).*exp(-1./(1-x(idx).^2)));

end

