function [ x,y,z,w ] = Gauss_LegendreSpherical(N)
% returns the set of quadrature points on which to evaluate the function
% and the set of weights to be use in the quadrature



thetas = (0:(2*N-1))*pi/N;
[x1D,w1D] = Gauss_Legendre1D(N,-1,1);
psis = acos(x1D);
w = w1D*2*pi;

thetas = thetas(:)';
psis = psis(:);

x = sin(psis)*cos(thetas);
y = sin(psis)*sin(thetas);
z = cos(psis)*(1+0*thetas);


end

