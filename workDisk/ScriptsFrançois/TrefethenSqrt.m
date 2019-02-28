function [ v ] = TrefethenSqrt(A,N,x,M,min_m,max_M)
%function [ v ] = TrefethenSqrt(A,u,N,x,M,min_m,max_M)
% Approximates v = sqrt(A)*u using trefethen's method.
% N is the number of iterations used.
if ~exist('M','var')||isempty(max_M)
    M = 1;
end

k2 = min_m/max_M; % elliptic functions parameter k^2
Kp = ellipke(1-k2);
t = 1i*(.5:N)*Kp/N;
[sn, cn, dn] = ellipj(imag(t),1-k2);
cn = 1./cn; dn = dn.*cn; sn = 1i*sn.*cn;
w = sqrt(min_m)*sn;
dzdt = cn.*dn;

if ~isempty(x)
    
    v = zeros(size(A,1),1);
    for j = 1:N
        v = v - dzdt(j)*((A-w(j)^2*M)\(M*x));
%        tmp = Sherman_MorissonInv(A-w(j)^2*M,u,M*x);
%        v = v - dzdt(j)*tmp;
    end
    v = (-2*Kp*sqrt(min_m)/(pi*N))*(M\(A*v));
%    v = (-2*Kp*sqrt(min_m)/(pi*N))*(M\(A*v + (u'*v)*u));
    v = M*v;
else
    
    Mat = zeros(size(A));
    for j = 1:N
        Mat = Mat - dzdt(j)*inv(A-w(j)^2*M)*M;
    end
    Mat = (-2*Kp*sqrt(min_m)/(pi*N))*inv(M)*(A*Mat);
    Mat = M*Mat;
    v = Mat;
end
end

function[v] = Sherman_MorissonInv(A,u,x)
% Computes v such that (A + u*u')v = x where u is rank 1

tmp1 = A\x;
tmp2 = A\u;
v = tmp1 - (tmp2*(u'*tmp1))/(1 + u'*tmp2);


end

