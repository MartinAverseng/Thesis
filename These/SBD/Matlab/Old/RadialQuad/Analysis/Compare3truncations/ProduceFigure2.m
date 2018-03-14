%% Parameters 

clear all;
close all;

step = 80;
Pmin = 10;
Pmax = 1000;
Nmax = 6;

D = 2;  
acut1 = 0.03;
acut2 = 0.03;
tolSchmidt = 0;
c = Inf;
b = 1;
G = @(jean_baptiste)(log(jean_baptiste));
Gprime = @(bernard_l_hermite)(1./bernard_l_hermite);
t = linspace(acut2,1,1000);

figConv = figureConventions();

%% Loop over P


Pvec = Pmin:step:Pmax;
pp = 0;
for P = Pvec
pp = pp+1;

%% Method 3 : Gram-Schmidt

a = acut1;
tol = tolSchmidt;
[ alpha3,bessZs3,~] = computeBesselCoeffH1_bis(a,b,G,Gprime,tol,false,P,true);
close all

%% Method 1 : exact expansion

a = 0;
tol = 0;
[ alpha1,bessZs1,~] = computeBesselCoeffH1_bis(a,b,G,Gprime,tol,false,P,true );
close all

%% Method 2 : polynomial expansion

for N = 1:Nmax
    ak = computeAk(N,acut2);
ext = @(Alfred)((polyExt(G,acut2,1,ak,zeros(length(ak),1),Alfred)));
der = @(Emile)((polyExtDer(Gprime,acut2,1,ak,zeros(length(ak)-1,1),Emile)   ));

a = 0;
tol = 0;
[ alpha2{N},bessZs1,~] = computeBesselCoeffH1_bis(a,b,ext,der,tol,false,P,true );
close all
end

bessZs = BesselZeros(P,3*pi/4);
fApprox1 = @(x)(besselApprox(bessZs,alpha1,0,x,1));
err1(pp) = max(abs(fApprox1(t) - G(t)));

for N = 1:Nmax
    fApprox2{N} = @(x)(besselApprox(bessZs,alpha2{N},0,x,1));
    err2(N,pp) = max(abs(fApprox2{N}(t) - G(t)));
end

fApprox3 = @(x)(besselApprox(bessZs,alpha3,0,x,1));
try
    err3(pp) = max(abs(fApprox3(t) - G(t)));
catch
    err3(pp) = NaN;
end


end

save('err1');
save('err2');
save('err3');

%% Process error for the second method
err2 = min(err2,[],1);


%% Create the figures
nameMethods = figConv.nameMethods;

figure
semilogy(Pvec,(err1),'Linewidth',2,'DisplayName',nameMethods{1})
hold on
semilogy(Pvec,(err2),'Linewidth',2,'DisplayName',nameMethods{2})
semilogy(Pvec,(err3),'Linewidth',2,'DisplayName',nameMethods{3})
grid on
% Show the machine precision
% plot(Pvec,log(sqrt(eps))*ones(size(Pvec)),'LineWidth',2,'LineStyle','--','HandleVisibility','off');
% txt = 'square root of machine precision';
yl = ylim;
yl = [yl(1), yl(2) + 50];
ylim(yl);
% text(length(Pvec)/3,log(sqrt(eps)) + abs(diff(yl))/20,txt,'Fontsize',24,'HandleVisibility','off');
xlabel('Number of components $P$','Interpreter','LaTex');
ylabel('Logarithmic error $\log(\varepsilon)$','Interpreter','LaTex');
set(gca,'Fontsize',24);
box on
legend show
legend boxoff
matlab2tikz([figConv.texFolderPath 'evolutionError3method.tex']); 


