close all;
clear all;

P0 = 5;
a = 0.05;

rho0 = besselJroots(0,P0);
A0 = gramMatrix(a,1,rho0);
b0 = laplaceCoeffs(a,1,rho0);

alpha0 = A0\b0;

x = linspace(0.005,1,200);

x0 = linspace(0.005,0.999,60);
vals0 = coeffTofunc(alpha0,rho0,x0);

plot(x,log(x),'k--','LineWidth',1,'DisplayName','$\log(r)$');
hold on

plot(x0,vals0,'-','DisplayName','SBD, $P=5$');
legend show;
xlabel('$r$');
box off;
axis tight;
legend('location','southeast');
legend boxoff;

set(gca,'XTick',[0.05 0.5 0.75 1]);
set(gca,'XTickLabel',{'$a=0.05$', '0.5', '0.75', '1'});
currentDir = fileparts(mfilename('fullpath'));
path = fullfile(currentDir,'LogVsSBD.tex');
matlab2tikz(path,'width','0.9\plotwidth','parseStrings',false,...
    'extraTikzpictureOptions','trim axis left, trim axis right');

