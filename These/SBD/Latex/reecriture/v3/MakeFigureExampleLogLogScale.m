close all;
clear all;

P0 = 5;
P1 = 30;
P2 = 50;
a = 0.05;

rho0 = besselJroots(0,P0);
A0 = gramMatrix(a,1,rho0);
b0 = laplaceCoeffs(a,1,rho0);

rho1 = besselJroots(0,P1);
A1 = gramMatrix(a,1,rho1);
b1 = laplaceCoeffs(a,1,rho1);

rho2 = besselJroots(0,P2);
A2 = gramMatrix(a,1,rho2);
b2 = laplaceCoeffs(a,1,rho2);

alpha0 = A0\b0;
alpha1 = A1\b1;
alpha2 = A2\b2;
x0 = linspace(0.005,0.999,60);
x1 = linspace(0.005,0.999,300);
x2 = linspace(0.005,0.999,500);
vals0 = abs(coeffTofunc(alpha1,rho0,x0)-log(x0));
vals1 = abs(coeffTofunc(alpha1,rho1,x1)-log(x1));
vals2 = abs(coeffTofunc(alpha2,rho2,x2)-log(x2));


semilogy(x0,vals0,'-o','LineWidth',0.7,'Markersize',4,'DisplayName','SBD, $P=5$');
hold on
semilogy(x1,vals1,'-x','LineWidth',0.7,'Markersize',4,'DisplayName','SBD, $P=30$');
semilogy(x2,vals2,'-s','LineWidth',0.7,'Markersize',4,'DisplayName','SBD, $P=50$');
hold on;
plot([a,a],ylim,'k--','HandleVisibility','off')
xlabel('$r$');
ylabel('$\abs{\log r - \sum_p\alpha_p e_p(r)}$')
box off;
axis tight;
set(gca,'XTick',[0.05 0.5 0.75 1]);
set(gca,'XTickLabel',{'$a=0.05$','0.5', '0.75', '1'});
set(gca,'YAxisLocation','right')
legend show;
legend('location','northoutside');
legend boxoff;

currentDir = fileparts(mfilename('fullpath'));
path = fullfile(currentDir,'LogVsSBDlogscale.tex');
matlab2tikz(path,'width','0.8\plotwidth','parseStrings',false,...
    'extraTikzpictureOptions','trim axis left, trim axis right');

