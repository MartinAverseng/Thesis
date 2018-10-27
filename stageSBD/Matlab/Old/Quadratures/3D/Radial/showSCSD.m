function[y] = showSCSD(Beta,Npoints,rho)
% Beta are the quadrature parameters for the SCSD method
P = length(Beta);
eps = exp(-2*sin(rho)*P);
r = linspace(0,pi,Npoints);
pp = 2*(0:P-1)'+1;
pp_r = pp*r;

bet = repmat(Beta,1,Npoints);
y = sum(bet.*sin(pp_r),1);

figure;
plot(r,y);
hold on;
title('SCSD approximation')
set(gca,'XTick',[rho, pi/2,pi - rho]);
set(gca,'XTickLabel',{'\rho', '\pi /2','\pi - \rho'});
plot(r,1+eps*ones(size(r)),'--r');
plot(r,1-eps*ones(size(r)),'--r');
plot([rho,rho],[0,1+eps],'--k');
plot([pi-rho,pi-rho],[0,1+eps],'--k');

figure;
plot(r,log(abs(1-y)));
title('Error for SCSD in dB')


end