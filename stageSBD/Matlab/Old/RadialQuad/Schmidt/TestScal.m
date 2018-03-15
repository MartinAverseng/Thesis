D = 3;
k = D/2-1;
a = sort(rand(1,2));
gamma = 1;

[u,w] = Gauss_Legendre1D(1000,a(1),a(2));

funcInt = @(x)(x.*besselj(k,gamma*x).^2);
if k==0
    J_kdiff = @(x)(-besselj(1,x));
else
    J_kdiff = @(x)(0.5*(besselj(k-1,x)-besselj(k+1,x)));
end
funcCroch = @(x)((x.^2-k^2).*besselj(k,x).^2 + x.^2.*J_kdiff(x).^2);

I1 = sum(w.*funcInt(u));
I2 = 1/(2*gamma^2)*(funcCroch(gamma*a(2)) - funcCroch(gamma*a(1)));