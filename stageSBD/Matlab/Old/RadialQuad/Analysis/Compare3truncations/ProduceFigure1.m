%% Parameters 

clear all;
close all;

D = 2;
acut = 0.2;
tolSchmidt = 1e-3;
c = Inf;
b = 1;
G = @(x)(log(x));
Gprime = @(x)(1./x);
figConv = figureConventions();

%% Method 3 : Gram-Schmidt

a = acut;
tol = tolSchmidt;
[ alpha3,bessZs3,~] = computeBesselCoeffH1_bis(a,b,G,Gprime,tol,false );
close all
Pmax = length(alpha3);

%% Method 1 : exact expansion

a = 0;
tol = 0;
[ alpha1,bessZs1,~] = computeBesselCoeffH1_bis( a,b,G,Gprime,tol,false,Pmax);
close all

%% Method 2 : polynomial expansion

N = 2;
ak = computeAk(N,acut);
ext = @(x)((polyExt(G,acut,1,ak,zeros(length(ak),1),x)));
der = @(x)((polyExt(Gprime,acut,1,(1:length(ak)-1)'.*ak(2:end),zeros(length(ak)-1,1),x)));


a = 0;
tol = 0;
[ alpha2,bessZs2,~] = computeBesselCoeffH1_bis( a,b,ext,der,tol,false,Pmax);

close all



%% Gather results


P = Pmax;
bessZs = BesselZeros(Inf,0,P,0);
fApprox1 = @(x)(besselApprox(bessZs,alpha1,0,x,1));
fApprox2 = @(x)(besselApprox(bessZs,alpha2,0,x,1));
fApprox3 = @(x)(besselApprox(bessZs,alpha3,0,x,1));

beta1 = alpha1(1:P)/norm(alpha1,2);
beta2 = alpha2(1:P)/norm(alpha2,2);
beta3 = alpha3(1:P)/norm(alpha3,2);


%% Create figures
close all;

nameMethods = figConv.nameMethods;

t = linspace(0,1,1000);
figure
plot(t,log(abs(fApprox1(t)-log(t))),'LineWidth',1,'DisplayName',nameMethods{1});
hold on
plot(t,log(abs(fApprox2(t)-log(t))),'LineWidth',1,'DisplayName',nameMethods{2});
hold on
plot(t,log(abs(fApprox3(t)-log(t))),'LineWidth',1,'DisplayName',nameMethods{3});
xlim([0,1])
yl = ylim;
ylim(yl);
xlabel('$|x|$','Interpreter','LaTex')
ylabel('Error (dB)')
set(gca,'Fontsize',24)
box on
legend show
legend boxoff 

matlab2tikz('C:\Users\Martin\Documents\Cours\SCSD\Latex\ComptesRendus\error3method.tex'); 

figure
scatter(bessZs,log(abs(beta1)),'LineWidth',2,'DisplayName','Truncated Fourier-Bessel of G')
hold on
scatter(bessZs,log(abs(beta2)),'LineWidth',2,'DisplayName','Truncated Fourier-Bessel of polynomial extension')
scatter(bessZs,log(abs(beta3)),'LineWidth',2,'DisplayName','Gram-Schmidt')
xlim([0,bessZs(end)+pi])
yl = ylim;
ylim([yl(1),yl(2)+4]);
xlabel('$\rho_p$','Interpreter','LaTex')
ylabel('$\log(|\alpha_p|)$','Interpreter','LaTex')
set(gca,'Fontsize',24)
box on
legend show
legend boxoff 
matlab2tikz('C:\Users\Martin\Documents\Cours\SCSD\Latex\ComptesRendus\spectrum3method.tex'); 

