%% Example du noyau G(r) = 1/r
clear all;
close all;
col = get(gca,'ColorOrder');
nmax= 5;
lambda = [1 1 9 9*25 9*25*49 9*25*49*81];
n = 0;
M = multiDirMatrix(n);
mu = M\lambda(1:n+1)';
P = 300;
gamma = linspace(0.5,7,5);
markers = {'-o','-x','-s','-*','-d'};

rho = besselJroots(0,P);
for j = 1:length(gamma)
    a = gamma(j)/P;
    A = gramMatrix(a,1,rho);
    b = invxCoeff(a,rho);
    b2 = b;
    for k = 1:n
        b2 = b2 - mu(k+1)*xpow2nCoeffs(a,rho,k);
    end
    alpha = A\b;
    alpha2 = A\b2;
    x = linspace(0,0.99,100);
    target = (1./(x));
    vals = coeffTofunc(alpha,rho,x)+1;
    vals2 = coeffTofunc(alpha2,rho,x);
    for k=0:n
        vals2 = vals2 + mu(k+1)*x.^(2*k);
    end
    semilogy(x,abs(target-vals),markers{j},'color',col(j,:),'DisplayName',sprintf('$\\Pa = %g$',gamma(j)));
    hold on;
    plot([a a],ylim,'--','color',col(j,:),'HandleVisibility','off');
end

xlabel('$r$');
ylabel('$Error magnitude$')
legend show;
legend boxoff;
box off;
axis tight;
set(gca,'XTick',[0.2, 0.4, 0.6, 0.8, 1]);

currentDir = fileparts(mfilename('fullpath'));
path = fullfile(currentDir,'arbitraryKernel1.tex');
matlab2tikz(path,'width','0.8\plotwidth','parseStrings',false,...
    'extraTikzpictureOptions','trim axis left, trim axis right');

%% Example du noyau G(r) = 1/r
clear all;
close all;
col = get(gca,'ColorOrder');
nmax= 5;
lambda = [1 1 9 9*25 9*25*49 9*25*49*81 9*25*49*81*121];
n = 3;
M = multiDirMatrix(n);
mu = M\lambda(1:n+1)';
P = 300;
gamma = linspace(0.5,7,5);
markers = {'-o','-x','-s','-*','-d'};

rho = besselJroots(0,P);
for j = 1:length(gamma)
    a = gamma(j)/P;
    A = gramMatrix(a,1,rho);
    b = invxCoeff(a,rho);
    b2 = b;
    for k = 1:n
        b2 = b2 - mu(k+1)*xpow2nCoeffs(a,rho,k);
    end
    alpha = A\b;
    alpha2 = A\b2;
    x = linspace(0,0.99,100);
    target = (1./(x));
    vals = coeffTofunc(alpha,rho,x)+1;
    vals2 = coeffTofunc(alpha2,rho,x);
    for k=0:n
        vals2 = vals2 + mu(k+1)*x.^(2*k);
    end
    semilogy(x,abs(target-vals2),markers{j},'color',col(j,:),'DisplayName',sprintf('$\\Pa = %g$',gamma(j)));
    hold on;
    plot([a a],ylim,'--','color',col(j,:),'HandleVisibility','off');
end

xlabel('$r$');
ylabel('$Error magnitude$')
legend show;
legend boxoff;
box off;
axis tight;
set(gca,'XTick',[0.2, 0.4, 0.6, 0.8, 1]);

currentDir = fileparts(mfilename('fullpath'));
path = fullfile(currentDir,'arbitraryKernel2.tex');
matlab2tikz(path,'width','0.8\plotwidth','parseStrings',false,...
    'extraTikzpictureOptions','trim axis left, trim axis right');

