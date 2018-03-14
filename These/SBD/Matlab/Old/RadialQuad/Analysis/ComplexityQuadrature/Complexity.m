%% Evaluate the complexity of the radial quadrature in N

%% Parametres : 
func = @(x)(log(x));
derivative = @(x)(1./x);
b = 1;
tol = exp(-5);
askGraph = false;

Ns = [500 800 1000 1500 2000 5000 8000 10000];
as = 1./Ns;
t = zeros(size(Ns));
P = zeros(size(Ns));
mem = zeros(size(Ns));
for iter = 1:length(Ns);
    fprintf('iteration number %d / %d \n',iter,length(Ns));
    tic;
    [~,rho,w0] = computeBesselCoeffH1_bis(as(iter),b,func,derivative,tol,askGraph);
    t(iter) = toc;
    P(iter) = length(rho);
    mem(iter) = monitor_memory_whos();
end

%% Figures : 
figConv = figureConventions;
% Number of components 
figure 
loglog(1./as,P,'o','HandleVisibility','off');
hold on
loglog(1./as,1./as+3.35,'--k');
grid on
h = legend({'$y = x + C$'},'Location', [0.7, 0.2, 0.11,0.3]);
set(h,'Interpreter','LaTex')
xlim([1/as(1),1./as(end)]);
set(gca,'FontSize',18);
ylabel('$\log(P)$','Interpreter','LaTex')
xlabel('$\log\left(\frac{1}{a}\right)$','Interpreter','LaTex')
box on
legend boxoff 
matlab2tikz([figConv.texFolderPath 'ComplexityNumberOfComponents.tex'])

% Computational time
figConv = figureConventions;
% Number of components 
figure 
loglog(1./as,t,'o','HandleVisibility','off');
hold on
loglog(1./as,100000*(1./as).^1*exp(-17),'--k');
grid on
h = legend({'$y = 3x + C$'},'Location', [0.7, 0.2, 0.11,0.3]);
set(h,'Interpreter','LaTex')
xlim([1/as(1),1./as(end)]);
set(gca,'FontSize',18);
ylabel('$\log(t)$','Interpreter','LaTex')
xlabel('$\log\left(\frac{1}{a}\right)$','Interpreter','LaTex')
box on
legend boxoff 
matlab2tikz([figConv.texFolderPath 'ComplexityTime.tex'])