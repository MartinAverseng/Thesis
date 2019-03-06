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
N   = 10
tol = 1e-4
typ = 'P1'

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
sigmaw = weightedDom(Disk2, 2);

% Finite elements on the sphere
Vh = fem(Disk2,typ);

% Incident wave
PW = @(X) ones(size(X,1),1);

% Mass matrix
tic
Mw = integral(sigmaw,Vh,Vh);
Mw = (Mw + Mw')/2;
tMw = toc
%[L,U,P,Q] = lu(Mw);
%invM = @(u)(Q*(U\(L \(P*u))));
invM = @(u)(Mw\u);



% Weighted Delta matrix 
weight = @(X) (1. - X(:,1).^2 - X(:,2).^2);
tic
Delta = integral(sigmaw,grad(Vh),weight,grad(Vh));
tDelta = toc


N = size(Mw,1);
lambda = 32/pi^3;
phi0 = ones(N,1);
diagM = Mw*phi0;
Mphi0 = sqrt(lambda)*diagM;

% Projected Green kernel  
Gxy = @(X,Y) femGreenKernel(X,Y,'[1/r]',0);

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' G(x,y) psi(y) dx dy 
tic
Sw = 1/(4*pi) .* integral(sigmaw,sigmaw,Vh,Gxy,Vh);
Srw  = 1/(4*pi) .* projRegularize(sigmaw, sigmaw, Vh, '[1/r]', Vh);
Sw = Sw + Srw;
tSw = toc

% Finite element incident wave trace --> \int_Sx psi(x)' pw(x) dx
RHSw = integral(sigmaw,Vh,PW);

%%
%tic
%sqr = TrefethenSqrt(Delta,5,[],Mw,1,3000);
%Tref = toc

Prec = @(u) (invM(TrefethenSqrt(Delta ,5,invM(u),Mw,0.7,1000)));
Prec2 = @(u) (invM(Delta*invM(Sw*u) + Mphi0*dot(Mphi0,invM(Sw*u))));

t1 = tic;
[lambda1,flag,relres,it,res1] = gmres(Sw,RHSw,100,1e-8,3000,Prec);
t1 = toc(t1);

t3 = tic;
[lambda1,flag,relres,it,res3] = gmres(Sw,RHSw,100,1e-8,3000,Prec2);
t3 = toc(t3);

t2 = tic;
[lambda2,flag,relres,it,res2] = gmres(Sw,RHSw,100,1e-8,3000);
t2 = toc(t2);

N1 = size(res1,1);
N2 = size(res2,1);
N3 = size(res3,1);
plot([1:N1]',log10(res1/norm(Prec(RHSw))),'g',[1:N2]',log10(res2),'r',[1:N3]',log10(res3/norm(Prec2(RHSw))),'k')
legend(num2str(t1),num2str(t2),num2str(t3));
% figure
% m = Disk2;
% m.vtx(:,3) = lambda1;
% plot(m,lambda1)
% view(0,90)
% axis equal;