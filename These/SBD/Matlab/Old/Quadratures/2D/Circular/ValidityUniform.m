%% Figure jusqu'à quand l'approx quad circulaire est valide


M = 30;
x = linspace(0,50,500);


figure
figConv = figureConventions;
plot(x,besselj(0,x),'LineWidth',2);
for i = 1:length(x)
    Japprox(i) = approxJ0circular(M,x(i),0);
end
hold on
plot(x,Japprox,'r--','LineWidth',2);
set(gca,'FontSize',15);
legend({'$J_0(x)$','Circular quadrature'},'Interpreter','LaTex');
legend boxoff
box on
grid on;
xlabel('$x$','Interpreter','LaTex')
matlab2tikz([figConv.texFolderPath '\circQuadUnifApprox.tex']);

figure
semilogy(x,(abs(besselj(0,x)-Japprox)),'LineWidth',2);
set(gca,'FontSize',15);
box on
grid on;
xlabel('$x$','Interpreter','LaTex');
ylabel('Quadrature error');
matlab2tikz([figConv.texFolderPath '\circQuadUnifApproxLogErr.tex']);