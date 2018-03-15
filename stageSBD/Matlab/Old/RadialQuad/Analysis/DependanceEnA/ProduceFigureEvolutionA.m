%% ProduceFigure Evolution A
clear all
close all
%% Parameters 
G = @(X)(log(X));
Gprime = @(X)(1./X);
b = 1;
askGraph = false;
aMin = 0.001;
aMax = 0.5;
Na = 50;
as = exp(linspace(log(aMax),log(aMin),Na));
tols = [1e-1 1e-3 1e-6];
i = 0;
figure
colorBase = [0.4660, 0.6740, 0.1880].^(1);

%% I°) Dépendance de P_{\varepsilon} en 1/a
for tol = tols
    j = 0;
    i = i+1;
    for a = as
        j = j+1;
        [~,rho,w0] = computeBesselCoeffH1_bis(a,b,G,Gprime,tol,askGraph);
        P(i,j) = length(rho);
    end
    hold on
    loglog(1./as,P(i,:),'color',colorBase.^(2*length(tols) +1  - 2*i),'LineWidth',2,'DisplayName',sprintf('tol = %s',num2str(tol)));
    hold off
end
xlabel('$\frac{1}{a}$','Interpreter','LaTex');
xlim([0 1350]);
grid on
ylabel('Number of components $P$','Interpreter','LaTex')
set(gca,'Fontsize',24);
box on
legend show
legend boxoff 
figConv = figureConventions;
matlab2tikz([figConv.texFolderPath 'evolutionOfA.tex']); 

%% II°) Evolution de k(\varepsilon)

clear all
close all

Ntol = 50;
tolMin = 1e-9;
tolMax = 1;
tols = exp(linspace(log(tolMax),log(tolMin),Ntol));

as = [0.008, 0.007,0.006, 0.005];
i = 0;
for tol = tols
    j = 0;
    i = i+1;
    for a = as
        j = j+1;
        [~,rho,w0] = computeBesselCoeffH1_bis(a,b,G,Gprime,tol,askGraph);
        P(i,j) = length(rho);
    end
    k(i) = max(P(i,:).*as);
end

figConv = figureConventions;
figure
plot(log(1./tols),k,'LineWidth',2,'HandleVisibility','off')
hold on
semilogx(log(1./tols),0.3*-log(tols)+0.14,'LineWidth',2,'LineStyle','--','color','k');
semilogx(log(1./tols),0.29*-log(tols)-0.6,'LineWidth',2,'LineStyle','--','color','k');
ylim([-2,8]);
xlim([0 23]);
legend({'$\gamma_1 = 0.3$, $\gamma_2 = 0.14$','$\gamma_1 = 0.29$, $\gamma_2 = -0.6$'},'Interpreter','LaTex');

grid on
legend boxoff
box on
set(gca,'FontSize',24);
xlabel('$-\log(\varepsilon)$','Interpreter','LaTex');
ylabel('Estimated $k(\varepsilon)$','Interpreter','LaTex');
matlab2tikz([figConv.texFolderPath '/upperLowerBoundk.tex']);