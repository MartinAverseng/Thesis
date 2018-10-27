%clear all
%close all
format SHORTENG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Résolution de
%  | Omega = [0,2pi]^2
%  | u - Laplacien(u) = f dans Omega
%  | du/dn = 0 sur le bord de Omega
%  | f(x,y) = cos(x)cos(y)
% solution exacte : u(x,y) = 1/3 cos(x)cos(y)
% -> errL2 = int((u-uh)^2)
% -> errH1 = int((u-uh)^2+(dx(u)-dx(uh))^2+(dy(u)-dy(uh))^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

type = 'Lagrange';
formuleInt = 5; % pour être exact en P4 (P2xP2)
errL2P0=[];
errL2P1=[];
errH1P2=[];
errH1P0=[];
errH1P1=[];
errL2P2=[];
h=[];

for nbPts=5:5:100
    tic
    %domaine [O,Lx]x[0,Ly]
    Lx = 2*pi;
    Ly = 2*pi;
    g.points=[ [(0:Lx/nbPts:Lx)' zeros(nbPts+1,1)] ; ...
               [Lx*ones(nbPts,1) (Ly/nbPts:Ly/nbPts:Ly)']; ...
               [((Lx-Lx/nbPts):(-Lx/nbPts):0)' Ly*ones(nbPts,1)];...
               [zeros(nbPts-1,1) ((Ly-Ly/nbPts):(-Ly/nbPts):Ly/nbPts)'] ];
    nbPoints = size(g.points,1);
    g.segments=[[(1:(nbPoints-1))' (2:nbPoints)'] ; [nbPoints 1]];
    
    aireMax=1.5*(( max(g.points(:,1))-min(g.points(:,1)) )/nbPts)^2;
    m=triangle(g,aireMax);
    
    %structures d'intégration
    %integ = buildInteg(m,formuleInt);
    %Nint = length(integ.w);
    %tic
    %W=spdiags(integ.w,0,Nint,Nint);
    %t0=toc
    % structures d'éléments finis

    efP1 = FEspace( m, type, 1, formuleInt);
    efP2 = FEspace( m, type, 2, formuleInt);
    % CALCUL DE h
    h = [h, mean(sqrt((m.vertices(m.edges(:,1),1)-m.vertices(m.edges(:,2),1)).^2 ...
              + (m.vertices(m.edges(:,1),2)-m.vertices(m.edges(:,2),2)).^2))];
    
    % Matrice du problème
    AP1 = efP1.u' * efP1.W * efP1.u + efP1.dxu' * efP1.W * efP1.dxu + efP1.dyu' * efP1.W * efP1.dyu;
    AP2 = efP2.u' * efP2.W * efP2.u + efP2.dxu' * efP2.W * efP2.dxu + efP2.dyu' * efP2.W * efP2.dyu;
    % Second membre
    fint = cos(efP1.integ.x).*cos(efP1.integ.y);
    FP1 = efP1.u' * efP1.W * fint;
    FP2 = efP2.u' * efP2.W * fint;

    % Résolution
    UhP1ddl = AP1\FP1;
    UhP2ddl = AP2\FP2;
    
    t=toc
    
    % Solution exacte aux pts d'intégration
    Uint   = cos(efP1.integ.x).*cos(efP1.integ.y)/3;
    dxUint = -sin(efP1.integ.x).*cos(efP1.integ.y)/3;
    dyUint = -cos(efP1.integ.x).*sin(efP1.integ.y)/3;
    
    % uh aux points d'intégration
    UhP1int = efP1.u * UhP1ddl;
    UhP2int = efP2.u * UhP2ddl;
    dxUhP1int = efP1.dxu * UhP1ddl;
    dxUhP2int = efP2.dxu * UhP2ddl;
    dyUhP1int = efP1.dyu * UhP1ddl;
    dyUhP2int = efP2.dyu * UhP2ddl;
    
    % erreur aux points d'intégration
    errintP1 = Uint-UhP1int;
    errintP2 = Uint-UhP2int;
    dxerrintP1 = dxUint-dxUhP1int;
    dxerrintP2 = dxUint-dxUhP2int;
    dyerrintP1 = dyUint-dyUhP1int;
    dyerrintP2 = dyUint-dyUhP2int;
    
    % erreur en norme L2 et H1
    errL2P1 = [errL2P1, sqrt(errintP1'*efP1.W*errintP1)];
    errL2P2 = [errL2P2, sqrt(errintP2'*efP2.W*errintP2)];
    errH1P1 = [errH1P1, sqrt(errintP1'*efP1.W*errintP1+dxerrintP1'*efP1.W*dxerrintP1+dyerrintP1'*efP1.W*dyerrintP1)];
    errH1P2 = [errH1P2, sqrt(errintP2'*efP2.W*errintP2+dxerrintP2'*efP2.W*dxerrintP2+dyerrintP2'*efP2.W*dyerrintP2)];
end


% figure(1)
% subplot(1,2,1)
% plotMesh(m)
% title('Solution P1, CB Neumann')
% hold on
% plotSol(m,efP1,integ,UhP1ddl,1)
% axis square
% hold off
% subplot(1,2,2)
% plotSol(m,efP1,integ,UhP1ddl,2)
% axis square
% hold off


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








