function [ A ] = d_gramMatrix( a,rho )
% Computes Dp Dq \int_{B(a,1)} D^2(e_p) : D^2(e_q)
% where e_p and e_q are the normalized eigenvectors of the laplacian on the
% unit ball for the H10 scalar product. 
% The constants Dp and Dq are re-normalization constants for the new scalar
% product, that is Dp = 1/||e_p||_2 with 
% ||e_p||_2 = \sqrt(\int_B(0,1) D^2(e_p):D^2(e_p))
rho = rho(:);
P = length(rho);
if length(a) == 2
    b = a(2);
    a = a(1);
else
    b = 1;
end
C = Cp(rho);
ek = @(r)(C.*besselj(0,rho*r));
el = @(r)(ek(r)');
dek = @(r)(-C.*rho.*besselj(1,rho*r));
del = @(r)(dek(r)');
rhok = rho;
rhol = rho';

B = gramMatrix([a,b],rho);
A = B*diag(rhol.^2) + 2*pi*a*dek(a)*(rhol.^2.*el(a)) + 2*pi*(dek(a)*del(a) - dek(b)*del(b));
%A = sqrt(diag(1./diag(A)))*A*sqrt(diag(1./diag(A)));

end

