clear all
close all
format SHORTENG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Résolution de
%  | Omega = [0,2pi]^2
%  | u - Laplacien(u) = f dans Omega
%  | u = 0 sur le bord de Omega
%  | f(x,y) = sin(x)sin(y)
% solution exacte : u(x,y) = 1/3 sin(x)sin(y)
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

for nbPts=3:3:50
    
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
    %plotMesh(m,'E');
    
    
    %structures d'intégration
    integ = buildInteg(m,formuleInt);
    Nint = length(integ.w);
    W=spdiags(integ.w,0,Nint,Nint);
    
    % structures d'éléments finis
    efP1 = FEspace( m, type, 1, formuleInt);
    efP2 = FEspace( m, type, 2, formuleInt);
    
    % CALCUL DE h
    h = [h, mean(sqrt((m.vertices(m.edges(:,1),1)-m.vertices(m.edges(:,2),1)).^2 ...
              + (m.vertices(m.edges(:,1),2)-m.vertices(m.edges(:,2),2)).^2))];
    
    % Gestion de la condition de Dirichlet via une matrice de prolongement
    % Matrice de prolongement P :
    %  - de taille Nddl x Nddl0
    %  - Si Uddl0 est définie sur V_h(Omega\bord) 
    %    alors Uddl = P * Uddl0 est définie sur V_h(Omega)
    %                           vaut Uddl0 sur l'intérieur de Omega
    %                           est nulle au bord de Omega
    % Construction de la matrice P pour les EF P1 :
    % - VI = indices des sommets du maillage internes à Omega
    % - Nddl0 = nb de ddl de V_h(Omega\bord) = length(VI)
    % - alors P(VI(j),j) = 1 pour j=1..Nddl0 et P est nulle ailleurs
    % Construction de VI :
    % - en P1 : VI = numéro comme ddl des sommets pour lesquels le label est le label interne
    %                (= numéro des sommets pour lesquels le label est le label interne)
    % - en P2 : VI = numéro comme ddl des sommets pour lesquels le label est le label interne
    %                (= numéro des sommets correspondants)
    %              + numéro comme ddl milieux des aretes dont un des deux sommets est interne 
    %                (= nb sommets + num arete correspondante, cf FEspace.m)
    VIP1 = find(m.lab_vertices == 0);
    Nddl0P1 = length(VIP1);
    NddlP1 = size(efP1.ddl,1);
    PP1 = sparse(VIP1',(1:Nddl0P1)',ones(Nddl0P1,1),NddlP1,Nddl0P1);
    VIP2 = [VIP1; find(m.lab_edges == 0) + size(m.vertices,1)];
    Nddl0P2 = length(VIP2);
    NddlP2 = size(efP2.ddl,1);
    PP2 = sparse(VIP2',(1:Nddl0P2)',ones(Nddl0P2,1),NddlP2,Nddl0P2);
    
    % Matrice du problème
    AP1 = efP1.u' * W * efP1.u + efP1.dxu' * W * efP1.dxu + + efP1.dyu' * W * efP1.dyu;
    A0P1 = PP1' * AP1 * PP1;
    AP2 = efP2.u' * W * efP2.u + efP2.dxu' * W * efP2.dxu + + efP2.dyu' * W * efP2.dyu;
    A0P2 = PP2' * AP2 * PP2;
        
    % Second membre
    fint = sin(integ.x).*sin(integ.y);
    FP1 = efP1.u' * W * fint;
    F0P1 = PP1' * FP1;
    FP2 = efP2.u' * W * fint;
    F0P2 = PP2' * FP2;
    
    % Résolution
    Uh0P1ddl = A0P1\F0P1;
    UhP1ddl = PP1 * Uh0P1ddl;
    Uh0P2ddl = A0P2\F0P2;
    UhP2ddl = PP2 * Uh0P2ddl;
        
    % Solution exacte aux pts d'intégration
    Uint   = sin(integ.x).*sin(integ.y)/3;
    dxUint = cos(integ.x).*sin(integ.y)/3;
    dyUint = sin(integ.x).*cos(integ.y)/3;
    
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
    errL2P1 = [errL2P1, sqrt(integ.w'*errintP1.^2)];
    errL2P2 = [errL2P2, sqrt(integ.w'*errintP2.^2)];
    errH1P1 = [errH1P1, sqrt(integ.w'*(errintP1.^2+dxerrintP1.^2+dyerrintP1.^2))];
    errH1P2 = [errH1P2, sqrt(integ.w'*(errintP2.^2+dxerrintP2.^2+dyerrintP2.^2))];

end

figure(1)
subplot(1,2,1)
plotMesh(m)
title('Solution P1, CB Dirichlet')
hold on
plotSol(m,efP1,integ,UhP1ddl,1)
axis square
hold off
subplot(1,2,2)
plotSol(m,efP1,integ,UhP1ddl,2)
axis square
hold off


figure(2)
subplot(1,2,1)
loglog(h,errL2P1,'b+-',h,errL2P2,'r+-',h,(h*10^-0.5).^2,'k--',h,(h*10^-0.5).^3,'k:')
legend('EF P1', 'EF P2','pente 2','pente 3')
title('erreur L2, CB Dirichlet')
xlabel('h')
subplot(1,2,2)
loglog(h,errH1P1,'b+-',h,errH1P2,'r+-',h,h,'k',h,(h*10^-0.5).^2,'k--')
legend('EF P1', 'EF P2','pente 1','pente 2')
title('erreur H1, CB Dirichlet')
xlabel('h')








