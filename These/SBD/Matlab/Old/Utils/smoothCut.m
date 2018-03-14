function [fc] = smoothCut(f,rho)

chi = @(r)(0.5+0.5*tanh(-1./(tan(pi*r)))); 
fc = @(x)((x<=rho).*chi(x/rho).*f(x) + and(x>=rho,x<=1-rho).*f(x) + (x>1-rho).*(1-(1-chi((1-x)/rho))).*f(x));

end

