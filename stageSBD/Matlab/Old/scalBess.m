function[v] = scalBess(lamj,lamk,a,b)


if lamj~=lamk
   t = @(x,lamk,lamj)(x*besselj(0,lamk*x)*lamj*besselj(1,lamj*x));
   u = @(x)(t(x,lamk,lamj) - t(x,lamj,lamk));
   v = (u(b) - u(a))/(lamj^2 - lamk^2);
else
   t = @(x)(x^2/2*(besselj(0,lamk*x)^2 + besselj(1,lamk*x)^2));
   v  = t(b) - t(a);
    
end


end