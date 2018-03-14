function [scal] = ScalBess(g1,g2,k,a,b)

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


end

