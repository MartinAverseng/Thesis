function Bessel
%clear all
close all
addpath('lib/');
tic
Lx = 2*pi;
N = 200;
theta = (0:2*pi/N:2*pi-2*pi/N)';

g.points=[cos(theta) sin(theta)];
g.segments=[(1:N-1)' (2:N)';N 1];
aireMax=0.2*(Lx/N)^2;
disp(['Geometrie ',num2str(toc),' s']);

tic
m = mesh2D(g,aireMax);
disp(['Maillage du domaine ',num2str(toc),' s']);
tic
%plotMesh(m,'E');
type = 'Lagrange';
efP1 = FEspace(m, type, 2);
[AP1,FP1] = BesselBil(m, efP1, efP1);
disp(['Temps d''assemblage ',num2str(toc),' s']);

tic
% Résolution
Uh1 = solve(AP1, FP1);
disp(['Temps de résolution ',num2str(toc),' s']);

% Affichage
plotSol(m,efP1,Uh1,1)
end

function [A,F] = BesselBil(m, test, trial)
A = intg(2, m, 'test', Id(test), 'trial', Id(trial)) ...
    + intg(2, m, 'test', Dx(test), 'trial', Dx(trial)) ...
    + intg(2, m, 'test', Dy(test), 'trial', Dy(trial)) ...
    + intg(1, m, 'label', 1, 'test', Id(test), 'trial', Id(trial));

un = @(x) ones(size(x,1),1);
F = intg(1, m, 'label', 1, 'test', Id(test), 'func', un);
end
