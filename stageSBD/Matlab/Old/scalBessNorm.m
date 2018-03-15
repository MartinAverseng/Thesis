function [ scal ] = scalBessNorm(g1,g2,a,b)
k = 0;

if k==0
    Jk_diff = @(x)(-besselj(1,x));
else
    Jk_diff = @(x)(0.5*(besselj(k-1,x)-besselj(k+1,x)));
end

if g1~=g2
    func1 = @(gg1,gg2,x)(gg2*x.*Jk_diff(gg2*x).*besselj(k,gg1*x));
    func2 = @(x)(func1(g1,g2,x)-func1(g2,g1,x));
    scal = 1/(g1^2-g2^2)*(func2(b)-func2(a));
else
    g = g1;
    funcCroch = @(x)((x.^2-k^2).*besselj(k,x).^2 + x.^2.*Jk_diff(x).^2);
    scal = 1/(2*g^2)*(funcCroch(g*b) - funcCroch(g*a));
end

scal = scal*2/abs(besselj(1,g1)*besselj(1,g2));
% Check
% [x,w] = Gauss_Legendre1D(4000,a,b);
% check = 2*sum(w.*x.*besselj(0,g1*x).*besselj(0,g2*x))...
%     /(abs(besselj(1,g1)*besselj(1,g2)));
% err = abs(scal - check);
% assert(err<1e-12,'Problem in computed scalar product');
end

% Validated
