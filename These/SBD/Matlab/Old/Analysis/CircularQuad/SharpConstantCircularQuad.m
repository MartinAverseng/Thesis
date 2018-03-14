%% Recherche d'une constante sharp pour l'approximation circulaire


clear all
close all
%% I°) Dépendance en M
rMin = 1;
rMax = 100000;
Nr = 100;
rs = linspace(rMin,rMax,Nr);
i = 0;
tol = 1e-12;
phis = 0:0.01:(2*pi);

for r = rs
    i = i+1
    reachedTol = false;
    M = fix(r + 12*r.^(0.3));
    while ~reachedTol
        M = M+1;
        reachedTol = true;
        for phi = phis
            x1 = r *cos(phi);
            x2 = r *sin(phi);
            errPhi = abs(approxJ0circular(M,x1,x2)-besselj(0,r));
            reachedTol = and(reachedTol,errPhi < tol);
            if ~reachedTol
                break
            end
        end
    end
    Ms(i) = M;
    
end

% Ms ~ rs + alpha x rs ^ beta
figure
loglog(rs,Ms./rs-1)
hold on

alpha = 13;
beta = 0.33;

figConv = figureConventions;
figure
loglog(rs,Ms-rs,'LineWidth',2);
hold on
loglog(rs,alpha*rs.^beta,'k--','LineWidth',2);
grid on;
h = legend({'$M_2(r)-r$','$\gamma_3 r^{\gamma_4}$'},'Interpreter','LaTex');
set(h,'Interpreter','LaTex');
set(gca,'FontSize',12);
xlabel('$r$','Interpreter','LaTex');
legend boxoff
box on
matlab2tikz([figConv.texFolderPath '/SharpCirc.tex'])