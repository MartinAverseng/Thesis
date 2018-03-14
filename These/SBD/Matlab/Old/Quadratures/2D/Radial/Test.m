%% Compute the development of a function in Fourier Bessel Series
close all;

Nr = 100;
rho = 0.1;
tol = 1e-6;
func = @(r)(-1/(2*pi)*log(r));
beta = besselQuad(func,rho,tol);
P = length(beta);
r01 = linspace(1/(5*Nr),1,Nr);
B = func(r01');

A = besselj(0,r01'*bessZs(1:P)');
figure

plot(r01,B,'DisplayName','Original function')
hold on; plot(r01,A*beta,'DisplayName','Approximation')
title(sprintf('Result of the quadrature, P = %d',length(beta)))
hold on; plot([rho rho],[0 B(1)],'k--','HandleVisibility','off');
set(gca,'XTick',[rho, 0.5 1]);
set(gca,'XTickLabel',{'\rho', 0.5, 1});
legend show;

figure 
err = (B - A*beta)';
plot(r01,log(abs(err)))
yl = ylim;
hold on; plot([rho rho],yl,'k--');
hold on; plot(r01,ones(size(r01))*log(tol),'r--');
title('Quadrature error (dB)');
set(gca,'XTick',[rho, 0.5 1]);
set(gca,'XTickLabel',{'\rho', 0.5, 1});

