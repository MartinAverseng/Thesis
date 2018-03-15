clear all;
close all;

%% Parameters

D = 2;
c = Inf;
a_vec = [0.02,0.05,0.1,0.2];
b = 1;
G = @(x)(log(x));
tol = -10i;
Pmax = 100;

%% Loop over $a$ and store the curves

for i = 1:length(a_vec)
    a = a_vec(i);
    [ ~,~,~,~,lambdaMin{i}] = BesselQuadSchmidt( D,c,G,a,b,tol,0,Pmax );
end
close all

%% Plot the curves

for i = 1:length(a_vec)
    figure
    plot(1./sqrt(lambdaMin{i}),'Linewidth',2);
    
    hold on;
    plot(pi*sqrt(1:100)/(1-a_vec(i)/2),'--','LineWidth',2);
    title(sprintf('a = %.2f',a_vec(i)));
    h = legend({'Numerical value','Conjecture'},'Location', [0.7, 0.3, 0.11,0.5]);
    set(h,'Interpreter','LaTex')
    legend show;
    legend boxoff
    
    set(gca,'Fontsize',24)
    box on
    
    path = fullfile('C:\Users\Martin\Documents\Cours\SCSD\Latex\ComptesRendus',sprintf('ConjnectureLambdaMin%d.tex',i));
    matlab2tikz(path);
end