% test de la majoration : 
clear all;
close all;

theta = 0.21;
r = 80;
z = r*exp(1i*theta);
x = real(z);
y = imag(z);
target = besselj(0,r);

for N = 1:200
    thet = (0:N-1)*2*pi/N;
    approx(N) = mean(exp(1i*real(z*exp(1i*thet))));
    err(N) = abs(approx(N) - target);
end

majoration1 = 2*(exp(1)*r./(1:N)).^(1:N);
figure
figConv = figureConventions;
semilogy(1:N,err,'LineWidth',2)
M1 = 80;
M2 = 140;
hold on
semilogy([M1 M1],([10^-16,2]),'k--','LineWidth',2);
semilogy([M2 M2],([10^-16,2]),'k--','LineWidth',2);
text(M1,10,'$M_1(r)$','Interpreter','LaTex','FontSize',15,'HorizontalAlignment','center');
text(M2,10,'$M_2(r)$','Interpreter','LaTex','FontSize',15,'HorizontalAlignment','center')
grid on;
xlabel('Number of points $M$','Interpreter','LaTex');
ylabel('Approximation error $e(r,M)$','Interpreter','LaTex')
set(gca,'FontSize',15)
ylim([1e-16,100])
matlab2tikz([figConv.texFolderPath '\errSuddenDrop_a.tex']);


