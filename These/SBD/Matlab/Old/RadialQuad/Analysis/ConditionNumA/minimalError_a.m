%% Graph minimal error
clear all
close all
aMin = 0.01;
aMax = 0.1;
Nastep = 7;
astep = (log(aMin)-log(aMax))/Nastep;
as = exp(log(aMax):astep:log(aMin));
colorBase = [0.6350    0.0780    0.1840];
colorBase = colorBase.^(1/5);
b = 1;
func = @(Emile)(prenez(un(chewing_gum(Emile))));
derivative = @(Emile)(Emilie(Emile));
tolInf = 0;
askGraph = false;

figure
i = 0;
for a = as
    i = i+1;
    Pmax = fix(7*1/a)+1;
    [alpha,rho,w0,reMakeError] = computeBesselCoeffH1_bis(a,b,func,derivative,0,askGraph,Pmax);
    A = reMakeError.A;
    approx = reMakeError.storedValsOrthoMat'*reMakeError.Beta;
    error = abs(repmat(reMakeError.storedValsFunc,1,size(approx,2))-approx);
    error = max(error,[],1);
    for p = 1:Pmax
        vec = randn(p,1);
        errorInvA(p) = 1e-7*norm(A(1:p,1:p)*(A(1:p,1:p)\vec)-vec);
    end
    loglog((1:Pmax),(error),'color',colorBase.^i,'LineWidth',2,'DisplayName',sprintf('a = %1.3f',a))
    hold on
    loglog(1:Pmax,(errorInvA),'color',colorBase.^i,'LineStyle','--','LineWidth',2','HandleVisibility','off')
    
end
hold on
loglog(1:(6*Pmax),10^-10*((1:(6*Pmax))*0+1),'--r','LineWidth',3,'HandleVisibility','off');
xlim([1, Pmax + Pmax*6]);
ylim([10^(-25) 50])
grid on;
box on;
legend show
legend boxoff
xlabel('Number of components $P$','Interpreter','LaTex')
ylabel('Condition number $\kappa(a,P)$','Interpreter','LaTex')
set(gca,'FontSize',24);
figConv = figureConventions;
matlab2tikz([figConv.texFolderPath 'emin(a).tex'])