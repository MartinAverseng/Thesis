%ScatteringDisk
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

clear all
close all

addpath('../Gypsilab/OpenMsh');
addpath('../Gypsilab/OpenMmg');
addpath('../Gypsilab/OpenDom');
addpath('../Gypsilab/OpenFem');
addpath('../Gypsilab/OpenHmx');

% Parameters
N   = 1e3
tol = 1e-4
typ = 'P1'
gss = 3;

% Meshes
Disk = mshDisk(N,1);
Disk.vtx(:,3) = sqrt(1.00000001 - Disk.vtx(:,1).^2 - Disk.vtx(:,2).^2);
[m,val] = mmg(Disk);
nrm = sqrt(sum(m.vtx.^2,2));
m.vtx(:,1) = m.vtx(:,1)./nrm;
m.vtx(:,2) = m.vtx(:,2)./nrm;
m.vtx(:,3) = m.vtx(:,3)./nrm;

Disk = m;
Disk.vtx(:,3) = 0;



% Domain
sigma = dom(Disk,gss);   

% Radiative mesh
square     = mshSquare(5*N,[5 5]);
square.vtx = [square.vtx(:,1) zeros(size(square.vtx,1),1) square.vtx(:,2)];
%plot(square)


% Incident wave
PW = @(X) ones(size(X,1),1);


%%% PREPARE OPERATOR
disp('~~~~~~~~~~~~~ PREPARE OPERATOR ~~~~~~~~~~~~~')

% Green kernel function --> G(x,y) = exp(ik|x-y|)/|x-y| 
Gxy = @(X,Y) femGreenKernel(X,Y,'[1/r]',0);

% Finite elements
Vh = fem(Disk,typ);

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' G(x,y) psi(y) dx dy 
tic
LHS = 1/(4*pi) .* integral(sigma,sigma,Vh,Gxy,Vh);
toc

% Regularization
tic
Sr  = 1/(4*pi) .* regularize(sigma,sigma,Vh,'[1/r]',Vh);
LHS = LHS + Sr;
toc

% Finite element incident wave trace --> \int_Sx psi(x)' pw(x) dx
RHS = integral(sigma,Vh,PW);

%%% SOLVE LINEAR PROBLEM
disp('~~~~~~~~~~~~~ SOLVE LINEAR PROBLEM ~~~~~~~~~~~~~')

% LU factorization
tic
[Lh,Uh] = lu(LHS);
toc

% Solve linear system [S] * lambda = P0
tic
lambda  = Uh \ (Lh \ RHS); % LHS \ RHS;
toc

% Graphical representation
figure
m = Disk;
m.vtx(:,3) = lambda.*sqrt(1.0000001 - m.vtx(:,1).^2 - m.vtx(:,2).^2);
plot(m,lambda)
axis equal;
%plot(square,abs(Pdom))
title('Total field solution')
colorbar


