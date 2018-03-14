function [y] = QuadVect(Beta,Npoints,a,b)
% Beta are the quadrature parameters for the SCSD method
% Return the approximation of unity using Beta model

P = length(Beta);
r = linspace(a,b,Npoints);
pp = 2*(0:P-1)'+1;
pp_r = pp*r;

bet = repmat(Beta,1,Npoints);
y = sum(bet.*sin(pp_r),1);

end

