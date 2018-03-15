%% Time RadialQuad
clear all
close all;
%% Parameters
a = 0.0001;
b = 1;
func = @(Emile)(log(Emile));
derivative = @(Emile)(1./Emile);
tolInf = 0;
checkCond = true;
askGraph = false;
Ps = fix(10*(1.5.^(1:20)))+1;
Pmin = 10;
Pmax = 8000;
NP = 50;
Ps = fix(exp(linspace(log(Pmin),log(Pmax),NP)));

i=0;
t = zeros(NP,1);
for P = Ps
    i = i+1;
    Pmax = P;
    
    t1 = tic;
    computeBesselCoeffH1_bis(a,b,func,derivative,tolInf,askGraph,Pmax,checkCond);
    t(i) = toc(t1);
end

%% Plot graph
figure
figConv = figureConventions;
loglog(Ps,t,'LineWidth',2,'HandleVisibility','off');
grid on
hold on
xlim([Ps(1) 10^4]);
ylim([10^(-5), 10^4]);
loglog(Ps(5:end),1/10000000*Ps(5:end).^2,'LineStyle','--','LineWidth',2);
legend({'$O(P^3)$'},'Interpreter','LaTex');
xlabel('Number of components $P$','Interpreter','LaTex');
ylabel('Computing time (s)');
legend boxoff
box on
set(gca,'FontSize',24);
matlab2tikz([figConv.texFolderPath 'ComplexQuadRad.tex']);