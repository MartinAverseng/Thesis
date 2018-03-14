close all;
clear all;

P0 = 5;
P1 = 30;
P2 = 50;
a = 0.05;

rho = besselJroots(0,P2);
A = gramMatrix(0.00001,1,rho);
b = laplaceCoeffs(0.00001,1,rho);
alpha = A\b;

rho0 = besselJroots(0,P0);
A0 = gramMatrix(a,1,rho0);
b0 = laplaceCoeffs(a,1,rho0);
alpha0 = A0\b0;

rho1 = besselJroots(0,P1);
A1 = gramMatrix(a,1,rho1);
b1 = laplaceCoeffs(a,1,rho1);
alpha1 = A1\b1;

rho2 = besselJroots(0,P2);
A2 = gramMatrix(a,1,rho2);
b2 = laplaceCoeffs(a,1,rho2);
alpha2 = A2\b2;

semilogy(rho,abs(alpha),'k-*','DisplayName','$c_p(G)$');
hold on;
semilogy(rho0,abs(alpha0),'-o','MarkerSize',4,'DisplayName','SBD : $P=5$');
semilogy(rho1,abs(alpha1),'Marker','x','DisplayName','SBD : $P=30$');
semilogy(rho2,abs(alpha2),'Marker','s','DisplayName','SBD : $P=50$');

legend show;
box off;
legend('location','southwest');
legend boxoff;
xlabel('$p$');
ylabel('$\alpha_p$');

currentDir = fileparts(mfilename('fullpath'));
path = fullfile(currentDir,'SBDCoeffLog.tex');
matlab2tikz(path,'width','0.9\plotwidth','parseStrings',false,'extraTikzpictureOptions',...
    'trim axis left, trim axis right'); 