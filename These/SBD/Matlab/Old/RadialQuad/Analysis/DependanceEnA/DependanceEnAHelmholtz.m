%% Parameters

clear all;
close all;

Pstep = 5;
Pmin = 5;
Pmax = 100;

a_step = 0.01;
a_min = 0.005;
a_max = 0.4;

D = 2;
tolSchmidt = 0;
c = Inf;
b = 1;
k = BesselYZeros(Inf,0,1,20);
G = @(x)(bessely(0,k*x));


figConv = figureConventions();

%% Loop over P and a


a_vec = a_min:a_step:a_max;
Pvec = Pmin:Pstep:Pmax;
errCompareAHelm = zeros(length(Pvec),length(a_vec));
pp = 0;
aa = 0;


for a = a_vec
    aa = aa+1;
    pp = 0;
    [ ~,bessZs,~,resL2,beta,B,normSQ ] = BesselQuadSchmidt( D,c,G,a,b,tolSchmidt,k,Pmax );
    
    for P = Pvec
        pp = pp+1;
        close all
        
        if length(beta) < P
            errCompareAHelm(pp,aa) = NaN;
        else
            alpha = B(1:P,1:P)*(normSQ(1:P,1:P)*beta(1:P));
            fApprox = @(x)(besselApprox(bessZs(1:P),alpha,0,x,0));
            t = linspace(a,1,1000);
            
            
            errCompareAHelm(pp,aa) = max(abs(fApprox(t) - G(t)));
        end
        
        
        
    end
end

save('errCompareAHelm');

%% Create the figures

figure
for anum = [1, 3, 5, 10, 16,34]
    a = a_vec(anum);
    hold on
    plot(Pvec,log(errCompareAHelm(:,anum)),'Linewidth',2,'DisplayName',sprintf('a = %.3f',a));
    % Rajouter un texte en bas de la courbe à côté
end
set(gca,'FontSize',24);
xlabel('Number of components P')
ylabel('Logarithmic error of approximation')
legend show
box on
legend boxoff
grid on
matlab2tikz([figConv.texFolderPath 'evolutionErrorInFunctionOfAHelmholtz.tex']); 



