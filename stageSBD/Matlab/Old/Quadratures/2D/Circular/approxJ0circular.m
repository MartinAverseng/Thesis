function [res] = approxJ0circular(M,x1,x2)

x1 = x1(:);
x2 = x2(:);

thetas = 2*pi/M*(0:M-1);
xis = exp(1i*thetas);
xi1 = real(xis);
xi2 = imag(xis);


res = sum(1/M*exp(1i*(x1*xi1 + x2*xi2)),2); 


end

