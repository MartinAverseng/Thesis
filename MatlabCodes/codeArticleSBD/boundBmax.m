function [ b ] = boundBmax(a,rho)

b = 2*pi*sum(F(rho*a)*2./(2*pi*rho.^2.*besselj(1,rho).^2));

    function[out] = F(u)
        out = u.^2/2.*(besselj(0,u).^2 + besselj(1,u).^2) - u.*besselj(0,u).*besselj(1,u);
    end


end

