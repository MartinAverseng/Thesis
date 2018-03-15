%% Estimation de D1 et D2
clear all;
close all;

Pvals = [50; 150; 500; 1500];
gamma = [linspace(0,1.470,10) linspace(1.6,10,25)];
markers = {'o','x','square','diamond'};
for i = 1:length(Pvals)
    P = Pvals(i);
    rho = besselJroots(0,P);
    for j = 1:length(gamma)
        a = gamma(j)/P;
        x = linspace(min(a*3,0.1),1,100);
        A = gramMatrix(a,1,rho);
        b = laplaceCoeffs(a,1,rho);
        alpha = A\b;
        testVals = linspace(a,max(min(5*a,1),0.1),100);
        quad = coeffTofunc(alpha,rho,testVals);
        Linf(i,j) = max(abs(log(testVals) - quad));
    end
    
    semilogy(gamma,Linf(i,:),markers{i},'MarkerSize',6,'DisplayName',sprintf('$P=%d$',P));
    hold on;
end
semilogy(gamma,5*exp(-3.7*gamma),'k-','DisplayName','$\Pa \mapsto 5\exp(-3.7\Pa)$');
hold on
semilogy(gamma,0*gamma + 10^(-10),'k--','HandleVisibility','off');
legend show
legend('location','best');
legend boxoff
box off
xlabel('$\gamma$');
ylabel('$L^{\infty}$ error');

set(gca,'FontSize',14);
set(gca,'YAxisLocation','right');
axis tight;
currentDir = fileparts(mfilename('fullpath'));
path = fullfile(currentDir,'ApplicationNumLaplace.tex');
matlab2tikz(path,'width','0.9\plotwidth','parseStrings',false,'addLabels',true,'extraTikzpictureOptions','trim axis left, trim axis right'); 