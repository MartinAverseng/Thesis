%% Laplace circle 

close all;
clear all;

%% Problem parameters :

tol = 1e-4;
r = 0.8;
Ncharges = 1000;
Nx = 500;
Ny = 500;
k = 23;

%% Geometry 

% Grid
Xmin = -1;
Xmax = 1;
Ymin = -1;
Ymax = 1;

Xaxis = linspace(Xmin,Xmax,Nx);
Yaxis = linspace(Xmin,Xmax,Ny)';
gridX = repmat(Xaxis,Ny,1);
gridY = repmat(Yaxis,1,Nx);
gridX = gridX(:);
gridY = gridY(:);
Ngrid = length(gridX);


% Charges
chargesX = r*cos(2*pi*1/Ncharges*(0:Ncharges-1));
chargesY = r*sin(2*pi*1/Ncharges*(0:Ncharges-1));
chargesQ = (rand(Ncharges,1)-1/2)*2/Ncharges;

% Plot the case geometry
figure
plot(chargesX,chargesY,'.');
xlim([Xmin,Xmax]);
ylim([Ymin,Ymax]);


%% Choose the parameter a
Rmax = 2*sqrt(Xmax^2 + Ymax^2);
bessYZs = BesselYZeros(Inf,0,2,k*Rmax);
Rmax = max(bessYZs)/k;
a = sqrt(sqrt(log(1/tol)))/sqrt(sqrt(Ncharges*Nx*Ny))/1.7; % To be modified later
Rmin = a*Rmax;

%% Compute the radial quadrature

D = 2;
c = Inf;
G = @(x)(bessely(0,k*x));
G_norm = @(X)(G(Rmax*X));
tolQuadRad = tol*a;

[ alpha,bessZs,~,resL2 ] = BesselQuadSchmidt( D,c,G_norm,a,1,tolQuadRad,k*Rmax,600);
fApprox = @(x)(besselApprox(bessZs,alpha,0,x,0));
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
sk = gridX/Rmax;
tk = gridY/Rmax;

field = nufft2d3(nj,xj,yj,cj,iflag,1e-4,nk,sk,tk)+ G(k*Rmax)*sum(chargesQ);

%% Produit proche

[I,rxy] = rangesearch([chargesX(:) chargesY(:)],[gridX(:) gridY(:)],1.5*Rmin);
jdx = cell2mat(I')';
rxy = cell2mat(rxy')';
idx = zeros(size(jdx));
j = 1;
for i=1:length(I);
    idx(j:j+length(I{i})-1) = i;
    j = j + length(I{i});
end

% Exact value of the green kernel
Gxy = G(rxy);
Gxy(rxy < 1e-15) = 0;
% SCSD values
SCSDxy = fApprox(rxy/Rmax) + G(Rmax);

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
set(gca,'Fontsize',24);
matlab2tikz('C:\Users\Martin\Documents\Cours\SCSD\Latex\ComptesRendus\CircleHelmholtz4.tex'); 


%% Comparaison avec la solution exacte en certains points
Ntest = 1000;
testInd = randsample(Ngrid,Ntest);

testField = zeros(Ntest,1);
for test = 1:Ntest
    t = testInd(test);
    for charge = 1:Ncharges
        r = sqrt((gridX(t)-chargesX(charge))^2+ (gridY(t)-chargesY(charge))^2 );
        testField(test) = testField(test) + bessely(0,k*r)*chargesQ(charge);
    end
end

err = norm(abs(field(testInd) - testField),'inf');
disp(err);

