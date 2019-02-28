%ScatteringDisk
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Assembler \omega \grad \grad + c\pi_0
% verifier les vecteurs propres harmoniques spheriques
% Mailler un disque irregulièrement.

clear all
close all

addpath('../Gypsilab/OpenMsh');
addpath('../Gypsilab/OpenMmg');
addpath('../Gypsilab/OpenDom');
addpath('../Gypsilab/OpenFem');
addpath('../Gypsilab/OpenHmx');

% Parameters
N   = 30
tol = 1e-4
typ = 'P1'
gss = 3;

% Mesh the disk and the half sphere S2
Disku = mshMyDisk(N,1);

% Ad hoc version

S2 = Disku;
r = sqrt(sum(S2.vtx.^2,2));
theta = pi/2*(1-r);
S2.vtx(2:end,1) = S2.vtx(2:end,1)./r(2:end).*cos(theta(2:end));
S2.vtx(2:end,2) = S2.vtx(2:end,2)./r(2:end).*cos(theta(2:end));
S2.vtx(:,3) = sin(theta(:));
Disk = S2;
Disk.vtx(:,3) = 0;
Disk2 = Disk;
Disk2.wgt = S2.ndv./Disk.ndv;

% Check weighted integral 
surface = sum(Disk2.ndv)

% Domain
sigmaw = dom(Disk2, gss);

% Finite elements on the sphere
Vh = fem(Disk,typ);

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' G(x,y) psi(y) dx dy 
weight = @(X) 1.000000001 - X(:,1).^2 - X(:,2).^2;
tic
Sw = integral(sigmaw,grad(Vh),weight,grad(Vh));
toc


% Mass matrix
tic
Mw = integral(sigmaw,Vh,Vh);
Mw = (Mw + Mw')/2;
toc

% Eigenvalues
[V,D,flag] = eigs(Sw, Mw, 16, 'smallestabs', 'IsSymmetricDefinite', true, 'Tolerance',1e-10);
%[V,D] = eig(full(Sw), full(Mw));
%[vp,I] =sort(diag(D));
% Graphical representation
figure
m = Disk;
for k = 1:16
    V(:,k) = V(:,k)/sqrt(V(:,k)'*Mw*V(:,k));
end
k = 0;
for i = 1:4
    for j = 1:4
        k = k+1;
        m.vtx(:,3) = V(:,k);
        subplot(4,4,k), plot(m,V(:,k))
        view(0,90)
        axis equal;
    end
end

