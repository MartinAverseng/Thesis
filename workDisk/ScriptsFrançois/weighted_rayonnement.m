%ScatteringDisk
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Assembler \omega \grad \grad + c\pi_0
% verifier les vecteurs propres harmoniques spheriques
% Mailler un disque irreguli�rement.

clear all
close all

addpath('../Gypsilab/OpenMsh');
addpath('../Gypsilab/OpenMmg');
addpath('../Gypsilab/OpenDom');
addpath('../Gypsilab/OpenFem');
addpath('../Gypsilab/OpenHmx');

% Parameters
N   = 20
tol = 1e-4
typ = 'P1'
gss = 3;

% % Incident wave
% PW = @(X) ones(size(X,1),1);

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
% Nrm = S2.nrm;
%Disk2.alpha = sum(S2.vtx(S2.elt(:,1),:).*Nrm,2);
% Check weighted integral 
surface = sum(Disk2.ndv)

% Domain
sigma  = dom(Disk,gss);   
sigmaw = dom(Disk2, gss);


% Projected Green kernel  
%Gxy = @(X,Y) proj1surR(X,Y);
Gxy = @(X,Y) femGreenKernel(X,Y,'[1/r]',0);

% Finite elements on the sphere
Vh = fem(Disk,typ);
x = [-1:0.001:1]';
theta = 0;
X = x*[cos(theta) sin(theta) 0];
X = Disk2.ctr;
X = X+ 0.0002*[ones(size(X,1),1),ones(size(X,1),1),zeros(size(X,1),1)];
% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' G(x,y) psi(y) dx dy 
tic
Sw = 1/(4*pi) .* integral(X,sigmaw,Gxy,Vh);
toc

% Regularization
tic
Srw  = 1/(4*pi) .* projRegularize(X, sigmaw, '[1/r]', Vh);
Sw = Sw + Srw;
toc
figure
%plot(x,Sw*ones(size(Sw,2),1))
plot(Sw*ones(size(Sw,2),1))
% % Finite element incident wave trace --> \int_Sx psi(x)' pw(x) dx
% RHSw = integral(sigmaw,Vh,PW);


