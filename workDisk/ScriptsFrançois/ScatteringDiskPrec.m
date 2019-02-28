%ScatteringDisk
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Assembler \omega \grad \grad + c\pi_0
% verifier les vecteurs propres harmoniques spheriques
% Mailler un disque irregulièrement.

%clear all
%close all

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

% Incident wave
PW = @(X) ones(size(X,1),1);

% Mesh the disk and the half sphere S2
Disku = mshMyDisk(N,1);
% dr = sqrt(pi/N);
% dr = 1/ceil(1/dr);
% dtheta = 2*pi/ceil(2*pi/dr); 
% Adapted mesh
% mmg version
%stp = Disku.stp;
%stp = dtheta;
%S2 = Disku;
%S2.vtx(:,3) = sqrt(1.00000000001 - Disku.vtx(:,1).^2 - Disku.vtx(:,2).^2);
%[S2,val] = mmg(S2,stp(3)*1.8);
%S2.vtx((S2.vtx(:,3)<1e-4),3) = 0;
%nrm = sqrt(sum(S2.vtx.^2,2));
%S2.vtx(:,1) = S2.vtx(:,1)./nrm;
%S2.vtx(:,2) = S2.vtx(:,2)./nrm;
%S2.vtx(:,3) = S2.vtx(:,3)./nrm;
%Disk = S2;
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
Nrm = S2.nrm;
Disk2.alpha = sum(S2.vtx(S2.elt(:,1),:).*Nrm,2);
%Disk2.alpha = ones(size(Disk2.alpha));
% Check weighted integral 
surface = sum(Disk2.ndv)

% Domain
sigma  = dom(Disk,gss);   
sigmau = dom(Disku,gss);   
sigmaw = dom(Disk2, gss);

% Radiative mesh
square     = mshSquare(5*N,[5 5]);
square.vtx = [square.vtx(:,1) zeros(size(square.vtx,1),1) square.vtx(:,2)];
%plot(square)

%%% PREPARE OPERATOR
disp('~~~~~~~~~~~~~ PREPARE OPERATOR ~~~~~~~~~~~~~')

% Projected Green kernel  
%Gxy = @(X,Y) proj1surR(X,Y);
Gxy = @(X,Y) femGreenKernel(X,Y,'[1/r]',0);

% Finite elements on the sphere
Vh = fem(Disk,typ);
Vhu = fem(Disku,typ);

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' G(x,y) psi(y) dx dy 
tic
Sw = 1/(4*pi) .* integral(sigmaw,sigmaw,Vh,Gxy,Vh);
S  = 1/(4*pi) .* integral(sigma,sigma,Vh,Gxy,Vh);
Su = 1/(4*pi) .* integral(sigmau,sigmau,Vhu,Gxy,Vhu);
toc

% Regularization
tic
Srw  = 1/(4*pi) .* projRegularize(sigmaw, sigmaw, Vh, '[1/r]', Vh);
Sr   = 1/(4*pi) .* regularize(sigma,  sigma,  Vh, '[1/r]', Vh);
Sru  = 1/(4*pi) .* regularize(sigmau, sigmau, Vhu,'[1/r]', Vhu);
Sw = Sw + Srw;
S  = S + Sr;
Su = Su + Sru;
toc

% Mass matrix
tic
Mw = integral(sigmaw,Vh,Vh);
M  = integral(sigma ,Vh,Vh);
Mu = integral(sigmau,Vhu,Vhu);
toc

% Finite element incident wave trace --> \int_Sx psi(x)' pw(x) dx
RHSw = integral(sigmaw,Vh,PW);
RHS  = integral(sigma,Vh,PW);
RHSu = integral(sigmau,Vhu,PW);

%%% SOLVE LINEAR PROBLEM
disp('~~~~~~~~~~~~~ SOLVE LINEAR PROBLEM ~~~~~~~~~~~~~')

%[lambda,flag,relres,it,res] = gmres(Sw,RHS,100,1e-6);
lambda1 = Sw\RHSw;
lambda2 = S\RHS;
lambda2 = lambda2.*sqrt(1.000000001 - Disk.vtx(:,1).^2 - Disk.vtx(:,2).^2);
lambdau = Su\RHSu;
lambdau = lambdau.*sqrt(1.000000001 - Disku.vtx(:,1).^2 - Disku.vtx(:,2).^2);
% Graphical representation
figure
m = Disk;
m.vtx(:,3) = lambda1;
subplot(1,3,1), plot(m)
view(0,20)
axis equal;

m.vtx(:,3) = lambda2;
subplot(1,3,2), plot(m)
view(0,20)
axis equal;

m = Disku;
m.vtx(:,3) = lambdau;
subplot(1,3,3), plot(m)
view(0,20)
axis equal;
%plot(square,abs(Pdom))
lambdaex = 4/pi*ones(size(lambda1));
lambdaexu = 4/pi*ones(size(lambdau));
error1 = sqrt((lambda1-lambdaex)'*M*(lambda1-lambdaex))/sqrt(lambdaex'*M*lambdaex)
error2 = sqrt((lambda2-lambdaex)'*M*(lambda2-lambdaex))/sqrt(lambdaex'*M*lambdaex)
erroru = sqrt((lambdau-lambdaexu)'*Mu*(lambdau-lambdaexu))/sqrt(lambdaexu'*Mu*lambdaexu)
size(lambda1)
size(lambdau)
%figure
%plot(log10(res))
