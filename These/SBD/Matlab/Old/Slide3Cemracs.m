%% Figure 2 

Ns = [1e3 5e3 1e4 5e4 1e5 5e5 1e6 5e6 1e7 5e7];
precomp = [0.045069 0.13875 0.26368 1.1734 2.3468 10.1695 22.8345 109.68 217.6512 1291.6096];
MV = [0.008488 0.038676 0.075192 0.36959 0.74023 4.0393 8.1473 46.45 87.9784 399.9699];

figure
loglog(Ns,precomp,'k-','LineWidth',2);
hold on
loglog(Ns,precomp,'ko','HandleVisibility','off')
loglog(Ns,MV,'k--','LineWidth',2)
loglog(Ns,MV,'ko','HandleVisibility','off')
hold on
loglog(Ns(7:end-1),Ns(7:end-1).*log(Ns(7:end-1))/300000,'b-','HandleVisibility','off','LineWidth',4);
legend({'pre-computation','MV product','$y \propto x \log(x)$'},'Interpreter','Latex');
xlabel('N','Interpreter','LaTex');
ylabel('t (s)','Interpreter','LaTex');
set(gca,'FontSize',18);
grid on
legend boxoff
axis tight