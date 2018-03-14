function Dirichlet
clear all
close all

addpath('lib/');

%domaine [O,Lx]x[0,Ly]
errL2P1=[]; errH1P2=[]; errH1P1=[]; errL2P2=[]; h=[];

%for nbPts=10:3:80
nbPts = 60;
    tic
    Lx = 2*pi;
    Ly = 2*pi;

    g.points=[0 0;Lx 0;Lx Ly;0 Ly];
    g.segments=[1 2;2 3;3 4;4 1];
    g.lab_segments=[3;3;2;2];
    aireMax=1.5*(Lx/nbPts)^2;
    
    m = mesh2D(g,aireMax);
    disp(['Maillage du domaine ',num2str(toc),' s']);
    plotMesh(m);

    sdfjsoidjf
    
    
    % Element fini P1 et P2
    type = 'Lagrange';
    efP1 = FEspace(m, type, 1);
    efP2 = FEspace(m, type, 2);
    
    [AP1, FP1, Dir1] = DirichletBil(m, efP1);
    [AP2, FP2, Dir2] = DirichletBil(m, efP2);
    
    % R�solution
    Uh1 = solve(AP1, FP1, Dir1);
    Uh2 = solve(AP2, FP2, Dir2);
    
    % Affichage
    % plotSol(m,efP1,Uh1,1)

    % CALCUL DE h
    h = [h, mean(sqrt((m.vertices(m.edges(:,1),1)-m.vertices(m.edges(:,2),1)).^2 ...
              + (m.vertices(m.edges(:,1),2)-m.vertices(m.edges(:,2),2)).^2))];
    
    % Solution exacte aux pts d'int�gration
    Uex = @(x) 0.1*x(:,1).*x(:,2) + sin(x(:,1)).*sin(x(:,2))/3;
    % erreur en norme L2 et H1
    errL2P1 = [errL2P1, erreur(m, Uex, Uh1, efP1, 'L2')];
    errL2P2 = [errL2P2, erreur(m, Uex, Uh2, efP2, 'L2')];
    errH1P1 = [errH1P1, erreur(m, Uex, Uh1, efP1, 'H1')];
    errH1P2 = [errH1P2, erreur(m, Uex, Uh2, efP2, 'H1')];
%end
sdjfsdj

figure(2)
subplot(1,2,1)
loglog(h,errL2P1,'b+-',h,errL2P2,'r+-',h,(h*10^-0.5).^2,'k--',h,(h*10^-0.5).^3,'k:')
legend('EF P1', 'EF P2','pente 2','pente 3')
title('erreur L2, CB Neumann')
xlabel('h')
subplot(1,2,2)
loglog(h,errH1P1,'b+-',h,errH1P2,'r+-',h,h,'k',h,(h*10^-0.5).^2,'k--')
legend('EF P1', 'EF P2','pente 1','pente 2')
title('erreur H1, CB Neumann')
xlabel('h')
end

function [A,F,D] = DirichletBil(m, ef)
A = intg(2, m, 'test', ef, 'trial', ef) ...
    + intg(2, m, 'test', Dx(ef), 'trial', Dx(ef)) ...
    + intg(2, m, 'test', Dy(ef), 'trial', Dy(ef));

f = @(x) 0.1*x(:,1).*x(:,2) + sin(x(:,1)).*sin(x(:,2));
F = intg(2, m, 'test', ef, 'func', f);
% Dirichlet condition
fdir = @(x) 0.1*x(:,1).*x(:,2);
D = DirCond(ef,[3,2],fdir);
end