close all;
clear all;

k = nextY0root(5);
P0 = 10;
a = 0.05;

rho0 = besselJroots(0,P0);
A0 = gramMatrix(a,1,rho0);
b0 = Y0Coeffs(a,1,rho0,k);

alpha0 = A0\b0;

x = linspace(0.005,1,500);

x0 = linspace(0.005,0.999,500);
vals0 = coeffTofunc(alpha0,rho0,x0);
hold on
plot(x,bessely(0,k*x),'k--','LineWidth',1,'DisplayName','$Y_0(kr)$, $k\approx 7.086...$');
hold on

plot(x0,vals0,'-','DisplayName','SBD, $P=10$');
legend show;
xlabel('$r$');
box off;
axis tight;
legend('location','southeast');
legend boxoff;

set(gca,'XTick',[0.05 0.5 0.75 1]);
set(gca,'XTickLabel',{'$a=0.05$', '0.5', '0.75', '1'});
currentDir = fileparts(mfilename('fullpath'));
path = fullfile(currentDir,'Y0VsSBD.tex');
matlab2tikz(path,'width','0.9\plotwidth','parseStrings',false,...
    'extraTikzpictureOptions','trim axis left, trim axis right');

