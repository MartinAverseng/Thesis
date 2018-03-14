%clear all
close all

addpath('lib/');

%domaine [O,Lx]x[0,Ly]

tic
Lx = 2*pi;
Ly = 2*pi;

g.points=[0 0;Lx 0;Lx Ly;0 Ly];
g.segments=[1 2;2 3;3 4;4 1];
g.lab_segments=[3;3;2;2];
aireMax=1.5*(Lx/nbPts)^2;

m = mesh2D(g,aireMax);
disp(['Maillage du domaine ',num2str(toc),' s']);
%plotMesh(m,'E');

% Physical constants
mu = 1;

% Element fini P1 et P2
type = 'Lagrange';
Mh = FEspace(m, type, 1);    % Pressure space
Vh = FEspace(m, type, 2, 2); % Velocity space of dim 2
ef = tensor(Vh,Mh);

% Integration
m.quad = 4; % pour être exact en P4 (P2xP2)

AP1 = mu*int2D(m, [], Dx(ef{1},1), Dx(ef{1},1)) ...
    + mu*int2D(m, [], (Dx(ef{1},2)+Dx(ef{2},1))/2, Dx(ef{1},2)) ...
    + mu*int2D(m, [], (Dx(ef{1},2)+Dx(ef{2},1))/2, Dx(ef{2},1)) ...
    + mu*int2D(m, [], Dx(ef{2},2), Dx(ef{2},2)) ...
    - int2D(m, [], ef{3}, Dx(ef{1},1)+Dx(ef{2},2)) ...
    - int2D(m, [], Dx(ef{1},1)+Dx(ef{2},2), ef{3}) ...

f = @(x,y) 0.1*x.*y + sin(x).*sin(y);
FP1 = int2D(m, f, Id(efP1));

% Dirichlet condition
fdir = @(x,y) 0.1*x.*y;
Dir1 = DirCond(efP1,[3,2],fdir);
% Résolution
Uh1 = solve(AP1, FP1, Dir1);


% Affichage
% plotSol(m,efP1,Uh1,1)


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
