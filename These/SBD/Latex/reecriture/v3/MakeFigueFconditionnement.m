close all;
clear all;
x2 = linspace(0,pi*1.471,40);
f = @(x2)(1-((1/2*(x2.^2+1)).*besselj(0, x2).^2-1/2*besselj(0, x2).*besselj(1, x2).*x2 + 1/2*besselj(1, x2).^2.*x2.^2-1/2));
estimate = f(x2);
plot(x2/pi,estimate,'k--','DisplayName','$F(\Pa)$');
markers = {'o','x'};
Pvals = [50; 500];
for i = 1:length(Pvals)
    P = Pvals(i);
    rho = besselJroots(0,P);
    for j = 1:length(x2)
        a2 = x2(j)/(pi*(P+1));
        A2 = gramMatrix(a2,1,rho);
        minEigA2(i,j) = min(eig(A2));
    end
    hold on
    plot(x2/pi,minEigA2(i,:),markers{i},'MarkerSize',6.5,'DisplayName',sprintf('$P=$%d',P));
end
legend show;
set(gca,'XTick',[0 0.5 1 1.471]);
set(gca,'YTick',[0 0.5 1]);
set(gca,'XTickLabel',{'0','0.5','1','$\Pastar$'});
xlabel('$\Pa$');
ylabel('$\lambda_{\min}(\gamma)$')

axis tight;
legend('location','southwest');
legend boxoff;
box off;
currentDir = fileparts(mfilename('fullpath'));
path = fullfile(currentDir,'Fconditionnement.tex');
matlab2tikz(path,'width','0.7\textwidth','height','5cm','parseStrings',false,'addLabels',true,'extraTikzpictureOptions','trim axis left, trim axis right'); 
