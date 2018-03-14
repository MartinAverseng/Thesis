close all;
clear all;

N = 1000;
x = 1:N;
rho = besselJroots(0,N);
C = 1./(sqrt(pi)*rho.*abs(besselj(1,rho)));
v = sqrt(2*pi*(1:N)').*C;
w = sqrt(2*pi*((1:N)'-1/4)).*C;

loglog(x,v-1,'o','MarkerSize',6,'DisplayName','$\sqrt{2\pi p}C_p - 1$');
hold on;
loglog(x,1-w,'x','MarkerSize',6,'DisplayName','$1 - \sqrt{2\pi (p-1/4)}C_p$');
legend show;

legend('location','southwest')
legend boxoff
xlabel('$p$');
set(gca,'FontSize',14);
currentDir = fileparts(mfilename('fullpath'));
path = fullfile(currentDir,'EncadremeentCp.tex');
box off; 
set(gca,'YTick',[1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1]);
set(gca,'XTick',[1,10,100,1000]);
matlab2tikz(path,'width','0.5\textwidth','height','5cm','parseStrings',false,'extraTikzpictureOptions','trim axis left, trim axis right'); 