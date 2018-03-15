%% ConditionNumber of A
clear all
close all

colorBase = [0.6350    0.0780    0.1840];
colorBase = colorBase.^(1/5);
b = 1;
Nstep = 15;
Pmin = 10;
Pmax = 1000;
step = (log(Pmax)-log(Pmin))/Nstep;
aMin = 0.001;
aMax = 0.2;
Nastep = 5;
astep = (log(aMin)-log(aMax))/Nastep;
as = exp(log(aMax):astep:log(aMin));
i = 0;
for a = as
    j = 0;
    i = i+1;
    A = matrixA(a,b,Pmax);
    for P = fix(exp(log(Pmin):step:log(Pmax)))
        j = j+1;
        B = A(1:P,1:P);
        condNum(i,j) = cond(B);
    end
    loglog( fix(exp(log(Pmin):step:log(Pmax))),condNum(i,:),'Color',colorBase.^i,'LineWidth',2,'DisplayName',sprintf('a = %1.3f',a));
    hold on;
end
xlim([Pmin, Pmax + Pmax*0.7*Pmin]);
grid on;
box on;
legend show
legend boxoff
xlabel('Number of components $P$','Interpreter','LaTex')
ylabel('Condition number $\kappa(a,P)$','Interpreter','LaTex')
set(gca,'FontSize',24);
figConv = figureConventions;
matlab2tikz([figConv.texFolderPath 'condNumA(a).tex'])
