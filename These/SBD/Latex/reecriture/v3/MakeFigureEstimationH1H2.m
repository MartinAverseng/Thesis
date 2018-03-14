% %% Test d'un truc
%
% clear all;
% close all;
% k = nextY0root(100);
% Pvals = [50; 150; 500; 1500];
% gamma = [linspace(0,1.470,10) linspace(1.6,10,25)];
%
% for i = 1:length(Pvals)
%     P = Pvals(i);
%     rho = besselJroots(k,P);
%     for j = 1:length(gamma)
%         a = gamma(j)/P;
%         A = gramMatrix(a,1,rho);
%         b = Y0Coeffs(a,1,rho,k) - 2/pi*.5772156649*J0Coeffs(a,1,rho,k);
%         alpha = A\b;
%         testVals = linspace(a,max(min(a+10*2*pi/k,1),0.1),500);
%         quad = coeffTofunc(alpha,rho,testVals);
% %         xplot = linspace(a,1,500);
% %         valplot = coeffTofunc(alpha,rho,xplot) + 2/pi*.5772156649*besselj(0,k*xplot);
% %         plot(xplot,valplot);
% %         hold on
% %         plot(xplot,bessely(0,xplot*k));
%         Linf(i,j) = max(abs(bessely(0,k*testVals) - 2/pi*.5772156649*besselj(0,k*testVals) - quad));
%     end
%     semilogy(gamma(2:end),Linf(i,2:end),'.','MarkerSize',6);
%     hold on;
% end




%% Graphe pour Helmholtz
clear all;
close all;

% 6 graphiques. La légende va à côté du graphique en haut à droite à
% l'extérieur.
ks = [10; 50; 100; 500; 1000; 5000];

markers = {'o','x','square','diamond'};

for kk = 1:6;
    close all;
    k = ks(kk);
    k = nextY0root(k);
    Pvals = [50; 150; 500; 1500];
    gamma = [linspace(0,1.470,10) linspace(1.6,10,25)];
    
    for i = 1:length(Pvals)
        P = Pvals(i);
        rho = besselJroots(k,P);
        for j = 1:length(gamma)
            a = gamma(j)/P;
            A = gramMatrix(a,1,rho);
            b = Y0Coeffs(a,1,rho,k);
            alpha = A\b;
            testVals = linspace(a,max(min(a+10*2*pi/k,1),0.1),500);
            quad = coeffTofunc(alpha,rho,testVals);
            %         xplot = linspace(a,1,100);
            %         valplot = coeffTofunc(alpha,rho,xplot);
            %         plot(xplot,valplot);
            %         hold on
            %         plot(xplot,bessely(0,xplot*k));
            Linf(i,j) = max(abs(bessely(0,k*testVals) - quad));
        end
        semilogy(gamma(2:end),Linf(i,2:end),markers{i},'MarkerSize',6,'DisplayName',sprintf('$P=%d$',P));
        hold on;
    end
    semilogy(gamma,5*exp(-3.7*gamma),'k-','DisplayName','$\Pa \mapsto 5\exp(-3.7\Pa)$');
    hold on
    semilogy(gamma,0*gamma + 10^(-10),'k--','HandleVisibility','off');
    set(gca,'XTick',[2,4,6,8]);
    if kk >=4
        xlabel('$\gamma$');
    else
        set(gca,'XTickLabels',[]);
    end
    if kk==1||kk==4
        ylabel('$L^{\infty}$ error');        
    else
        ylabel('');
        set(gca,'YTickLabels',[]);
    end
    title(sprintf('k=%.3f...',round(1000*k)/1000));
    currentDir = fileparts(mfilename('fullpath'));
    path = fullfile(currentDir,sprintf('HelmholtzConvSBD%d.tex',kk));
    axis tight;
    if kk==1
        legend show;
        matlab2tikz(path,'width','\linewidth','height','5cm','parseStrings',false,'extraTikzpictureOptions','trim axis left, trim axis right','extraAxisOptions','legend to name = legendref');
    else
        matlab2tikz(path,'width','\linewidth','height','5cm','parseStrings',false,'addLabels',false,'extraTikzpictureOptions','trim axis left, trim axis right');
    end
end




