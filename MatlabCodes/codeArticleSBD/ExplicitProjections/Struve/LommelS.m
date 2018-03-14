function [ out ] = LommelS( a,b,z)

out = z.^(a+1).*hypergeom(1, [3/2-(1/2)*b+(1/2)*a, 3/2+(1/2)*b+(1/2)*a], -(1/4)*z.^2)/((a-b+1)*(a+b+1));


end

