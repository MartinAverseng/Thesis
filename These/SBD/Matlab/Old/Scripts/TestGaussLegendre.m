%% Test the gauss-Legendre script

x1 = 1;
x2 = 0.3;
x3 = 0.01;
func = @(ksi1,ksi2,ksi3)(exp(1i*(x1*ksi1 + x2*ksi2 + x3*ksi3)));

Ns = 1:20;
err = zeros(size(Ns));
for N = Ns

I = integrateGauss(func,N);
target = 4*pi*sinc(sqrt(x1^2 + x2^2 + x3^2)/pi);

err(N) = abs(I - target);

end

plot(log(err));
