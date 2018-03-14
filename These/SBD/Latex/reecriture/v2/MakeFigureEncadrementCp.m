close all;
clear all;

N = 1000;
x = 1:N;
rho = besselJroots(0,N);
C = 1./(sqrt(pi)*rho.*abs(besselj(1,rho)));
v = sqrt(2*pi*(1:N)').*C;
w = sqrt(2*pi*((1:N)'-1/4)).*C;

semilogx(x,v-1,'.','MarkerSize',7);
hold on
semilogx(x,w-1,'.','MarkerSize',7);
semilogx(x,x*0,'k--','LineWidth',2);
xlabel('$p$');

set(gca,'FontSize',14);
set(gca,'XTick',[10.^(1:fix(log(N)/log(10))),N])
axis tight;
currentDir = fileparts(mfilename('fullpath'));
path = fullfile(currentDir,'EncadremeentCp.tex');
matlab2tikz(path,'width','0.5\textwidth','height','5cm','parseStrings',false,'addLabels',true,'extraTikzpictureOptions','trim axis left, trim axis right'); 
