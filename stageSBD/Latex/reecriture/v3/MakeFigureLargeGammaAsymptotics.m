clear all
close all;

gamma = [linspace(0,1.470,10) linspace(1.6,10,25)];
x1 = pi*gamma;
%estimate = 1 - (1/64)*x2.*(-x2.^3.*hypergeom([3/2, 2, 2], [1, 3, 3, 3], -x2.^2)+16*besselj(0, x2).^2.*x2+16*besselj(1, x2).^2.*x2-32*besselj(0, x2).*besselj(1, x2));
estimate = 1 - ((1/2*(x1.^2+1)).*besselj(0, x1).^2-1/2*besselj(0, x1).*besselj(1, x1).*x1 + 1/2*besselj(1, x1).^2.*x1.^2-1/2);
Pvals = [50; 150; 500; 1500];
h2 = figure;
semilogy(gamma,180*exp(-pi*1.85*gamma),'r','DisplayName','\autoref{conjLambdaMin}');
hold on;
semilogy(gamma,estimate,'k','DisplayName','\autoref{The:lowBoundCon}');
currentDir = fileparts(mfilename('fullpath'));
markers = {'o','x','square','diamond'};
try A = load(fullfile(currentDir,'minEigA1'));
    minEigA1 = A.minEigA1;
    x1temp = x1(1:end);
    for i = 1:length(Pvals)
        P = Pvals(i);
        gamma = x1temp/pi;
        semilogy(gamma,minEigA1(i,1:end),markers{i},'MarkerSize',5,'DisplayName',sprintf('$P=%d$',P));
    end
catch ex
    clear minEigA1;
    for i = 1:length(Pvals)
        P = Pvals(i);
        rho = besselJroots(0,P);
        for j = 1:length(x1)
            a1 = gamma(j)/P;
            A1 = gramMatrix(a1,1,rho);
            minEigA1(i,j) = min(eig(A1));
        end
        figure(h2)
        semilogy(gamma,minEigA1(i,:),markers{i},'MarkerSize',6,'DisplayName',sprintf('$P=%d$',P));
    end
    currentDir = fileparts(mfilename('fullpath'));
    save(fullfile(currentDir,'minEigA1'),'minEigA1');
end

legend show
legend('location','northeast');
legend boxoff;
hold on;
semilogy(gamma,0*gamma + eps,'k--');
box off;
[Xtick,I] = sort([1.471, get(gca, 'XTick')]); 
XtickLabel = get(gca, 'XTickLabel');
XticklabelNonSorted = {'$\Pastar$',XtickLabel{:}};
XTickLabel = XticklabelNonSorted(I);
set(gca, 'XTick',Xtick,'XTickLabel',XTickLabel);
set(gca,'YTick',10.^(-20:4:2));


xlabel('$\Pa$');
ylabel('$\lambda_{\min}(\gamma)$');
set(gca,'FontSize',14);
axis tight;
currentDir = fileparts(mfilename('fullpath'));
path = fullfile(currentDir,'minEigA1.tex');
matlab2tikz(path,'width','0.8\linewidth','height','5cm','parseStrings',false,'addLabels',true,...
    'extraTikzpictureOptions','trim axis left, trim axis right');
