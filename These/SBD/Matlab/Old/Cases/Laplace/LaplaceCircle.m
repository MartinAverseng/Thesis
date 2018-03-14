%% Laplace circle 

close all;
clear all;

%% Problem parameters :

tol = 1e-4;
r = 0.8;
Ncharges = 1000;


%% Geometry 

% Grid
Xmin = -1;
Xmax = 1;
Ymin = -1;
Ymax = 1;
Nx = 500;
Ny = 500;
gridX = linspace(Xmin,Xmax,Nx)';
gridY = linspace(Xmin,Xmax,Ny)';

% Charges
chargesX = r*cos(2*pi*1/Ncharges*(0:Ncharges-1));
chargesY = r*sin(2*pi*1/Ncharges*(0:Ncharges-1));
chargesQ = randn(Ncharges,1)*1/Ncharges;

% Plot the case geometry
figure
plot(chargesX,chargesY,'.');
xlim([Xmin,Xmax]);
ylim([Ymin,Ymax]);


%% Choose the parameter a
Rmax = 2*sqrt(2);
a = sqrt(sqrt(log(1/tol)))/sqrt(sqrt(Ncharges*Nx*Ny))/2; % To be modified later
Rmin = a*Rmax;

%% Compute the radial quadrature

D = 2;
c = Inf;
tt = 0;
b = 1-tt*a;
G = @(x)(log(x));
tolQuadRad = tol*sqrt(a);

[ alpha,bessZs,fApprox,resL2 ] = BesselQuadSchmidt( D,c,G,a,b,tolQuadRad );
P = length(alpha);

%% Circular quadrature 

asympt = log(8*norm(alpha,1)/tol);
ksi = [];
w = [];
for p = 1:P
    Np = fix((exp(1)*bessZs(p) + asympt))+1;
    err(p) = 3*(exp(1)*bessZs(p)/Np)^Np;
    ksi = [ksi bessZs(p)*exp(1i*(0:Np-1)/Np*2*pi)];
    w = [w alpha(p)/Np*ones(1,Np)];
end

%% Produit matrice vecteur lointain : 

tic

% 1. Passage en fourier non uniforme

nj = Ncharges;
xj = chargesX/Rmax;
yj = chargesY/Rmax;
cj = chargesQ;
iflag = -1; 

nk = length(ksi);
sk = real(ksi)';
tk = imag(ksi)';

aux = nufft2d3(nj,xj,yj,cj,iflag,1e-15,nk,sk,tk);

% 2. Retour en espace
nj = nk;
xj = sk;
yj = tk;
cj = aux.*w';
iflag = 1; 

nk = Nx*Ny;
sk = repmat(gridX',Ny,1)/Rmax;
tk = repmat(gridY,1,Nx)/Rmax;
sk = sk(:);
tk = tk(:);

field = nufft2d3(nj,xj,yj,cj,iflag,1e-15,nk,sk,tk);

%% Produit proche

[I,rxy] = rangesearch([chargesX(:) chargesY(:)],[sk(:)*Rmax tk(:)*Rmax],1.2*Rmin);
jdx = cell2mat(I')';
rxy = cell2mat(rxy')';
idx = zeros(size(jdx));
j = 1;
for i=1:length(I);
    idx(j:j+length(I{i})-1) = i;
    j = j + length(I{i});
end

% Exact value of the green kernel
Gxy = G(rxy/Rmax);
Gxy(rxy < 1e-15) = 0;
% SCSD values
SCSDxy = fApprox(rxy/Rmax);

% Matrix value
val = Gxy - SCSDxy;

% Compute correction matrix
M = sparse(idx,jdx,val,Nx*Ny,Ncharges);

field = field + M*chargesQ;


%% Affichage 

figure
ff = reshape(field,Nx,Ny);
imagesc(gridX,gridY,abs(ff));
xlabel('x')
ylabel('y')
matlab2tikz('C:\Users\Martin\Documents\Cours\SCSD\Latex\ComptesRendus\CircleLaplace.tex'); 


