% TinySCSD en 2 dimensions 
clear all;
close all;

tic;
% Création de deux nuages de points 
Nx = 100; % Number of charges
Ny = 100; % Nombre de points auxquels on souhaite connaître le champ
tol = 0.5; % Tolérance erreur. 

x = randn(Nx,2); % Vecteur position des charges, dans R^2
f = randn(Nx,1); % Valeur des charges (masse)
y = randn(Ny,2); % Vecteur position des points ou l'ont veut connaître le champ

% Affichage de la géométrie
figure
plot(z1,'*','DisplayName','Sources');
hold on;
plot(z2,'*','DisplayName','Targets');
title(sprintf('Points of the problem. Estimated radius : D = %.2f',D))
legend show

%% Définition de l'intéraction radiale G 
% On veut calculer qk = \sum_{l = 1..Nx}G(|y_k - x_l|)f_l



% Gid = 'log' % à utiliser pour résoudre une équation de Laplace
% Gid = 'Y0' % à utiliser pour résoudre une équation de Helmholtz
Gid = 'custom'; % Définition d'une fonction générale (nécessite de prendre de la marge à droite)
G = @(x)(1./x);


%% Définition de l'anneau d'approximation 

z1 = x(:,1) + 1i*x(:,2);
z2 = y(:,1) + 1i*y(:,2);

D = norm(abs(z1 - mean(z1)),'inf') + abs(mean(z1) - mean(z2)) + norm(abs(mean(z2)-z2),'inf'); % Taille du problème

% Choix de Rmin 
% Inclure ici la dichotomie. 
%Rmin = dichotomy();

% En attendant
Rmin = D/20;

% Choix de Rmax
switch(Gid)
    case 'log'
        Rmax = D;
    case 'Y0'
        Rmax = D;
    case 'custom'
        Rmax = D+Rmin; % On ajoute la même marge à la fin
        % On travaille sur l'anneau [Rmin, Rmax-Rmin]
    otherwise
        error('Kernel no supported')
end


%% Quadrature radiale

% Choix de la tolérace
% tolH1 = upperBound(tol); % A implémenter

% En attendant
delta = (Rmax - Rmin)/Nx; % petite cuisine supposant beaucoup de données 
% (proximité entre la norme L^2 et la moyenne des carrés aux points considérés. 

% Echelle normalisée : 
rho = Rmin/Rmax;
tolQuadRad = sqrt(delta*rho/(D*norm(q,2)))*tol/2;

% Fonction normalisée
G_norm = @(X)(G(Rmax*X) - G(Rmax));
% De sorte que G_norm(1) = 0, et G_norm(rho) = G(Rmin)
    
% Calcul de la décomposition de Fourier-Bessel
tic;
[beta,r_norm,quad,res] = besselQuadSchmidt(G_norm,rho,1,tolQuadRad);
% space = linspace(rho,1-rho,Nx+Ny);
% spaceErr = linspace(rho,1-rho,2*(Nx+Ny));
% reachedTol = false;
% p = 0;
% while ~reachedTol
%     p = p+1;
%     
%     zs = besselZeros(p);
%     RZ = zs*space;
%     A = besselj(0,RZ');
%     rhs = G_norm(space)';
%     weights = A\rhs;
%     err = norm(rhs - A*weights,1);
%     reachedTol = err < tol;
% end

tFourierBessel = toc;
P = length(beta); % number of components in quadrature

fprintf('Radial quadrature : %d components, \n %s seconds \n %d L2 error  \n\n',P,num2str(tQuad1),sqrt(res));


%% Quadrature circulaire 

tic

% Evaluation de l'erreur : 
asympt = log(8*norm(beta,1)*norm(q,'inf')/tol);

ksi = [];
w = [];
for p = 1:P
    Np = fix((exp(1)*r_norm(p) + asympt))+1;
    err(p) = 3*(exp(1)*r_norm(p)/Np)^Np;
    ksi = [ksi r_norm(p)*exp(1i*(0:Np-1)/Np*2*pi)];
    w = [w beta(p)/Np*ones(1,Np)];
end

tQuad2 = toc;
t1 = linspace(0,rho,1000);
func = raccordCn(rho,5);
hold on
plot(t1,func(t1)/Rmax-G(Rmax))
fprintf('Total quadrature : %d quadrature points, \n %s seconds \n\n',length(w),num2str(tQuad2));


%% Produit matrice vecteur lointain : 

tic

% 1. Passage en fourier non uniforme

nj = Nx;
xj = x(:,1)/Rmax;
yj = x(:,2)/Rmax;
cj = q;
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

nk = Ny;
sk = y(:,1)/Rmax;
tk = y(:,2)/Rmax;

field = nufft2d3(nj,xj,yj,cj,iflag,1e-15,nk,sk,tk) + G(Rmax)*sum(q);

tFFT = toc;
fprintf('SCSD product : \n %s seconds \n\n',num2str(tFFT));


%% Correction champ proche
tic

[I,rxy] = rangesearch(x,y,1.5*Rmin);
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
SCSDxy = quad(rxy/Rmax) + G(Rmax);

% Matrix value
val = Gxy - SCSDxy;

% Compute correction matrix
M = sparse(idx,jdx,val,Ny,Nx);
tSparse1 = toc;
disp(['Assemblage de la matrice de champ proche (s) : ',num2str(tSparse1)])

tic
field = field + M*q;
tSparse2 = toc;

disp(['Produit matrice-vecteur champ proche (s) : ',num2str(tSparse2)])


tSCSD = tQuad1 + tQuad2 + tSparse1 + tSparse2 + tFFT;

disp(['Total SCSD (s) : ',num2str(tSCSD)])



%% Comparaison solution exacte : 
% 
tic

trueField = zeros(Ny,1);
for iy = 1:Ny
    for ix = 1:Nx
        if sqrt((x(ix,1)-y(iy,1))^2 + (x(ix,2) - y(iy,2))^2)<1e-15
        else
            trueField(iy) = trueField(iy) + G(sqrt((x(ix,1)-y(iy,1))^2 + (x(ix,2) - y(iy,2))^2))*q(ix);
        end
        
    end
end

tExa = toc;

disp(['Solution exacte (s) : ',num2str(tExa)])





%% Error


error = norm(abs(field - trueField),'inf');
fprintf('Maximal error : %d \n',error)
