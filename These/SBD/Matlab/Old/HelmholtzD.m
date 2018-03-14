%-------------------------------------------------------------------------%
%      MyBEM 2 - Matthieu Aussal & Francois Alouges - Copyright 2016      %
%                                                                         %
% Ce logiciel MyBEM est la propriete de l'Ecole Polytechnique, tous les   %
% droits lui sont reserves, toute utilisation de ce logiciel est soumise  %
% a l'accord prealable et ecrit de l'Ecole Polytechnique.                 %
%                                                                         %
% This software MyBEM is owned by Ecole polytechnique, all rights are     %
% reserved, any use of this software is subjected to the written, prior   %
% consent of Ecole polytechnique.                                         %
%                                                                         %
% Contact :                                                               %
% matthieu.aussal@polytechnique.edu                                       %
% francois.alouges@polytechnique.edu                                      %
% martin.averseng@polytechnique.edu                                       %
%-------------------------------------------------------------------------%
%
% Creation : 2016.01.01
%
% Last Modification :
%
% Synopsis :

% Nettoyage ecran
clear all
close all
clc

% Initialisation
addpath('libGgNufft')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETRES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dimension du probleme
N   = 1e3;
tol = 1e-1;
k   = 15;

% Nuages de points dans un cube
X = -1 + 2*rand(N,2);
Y = -1 + 2*rand(N,2);

% % Nuages de points sur un cerc
% theta = (1:N)' .* (2*pi/N);
% X = [cos(theta),sin(theta)];
% X = unique(X,'rows');
% Y = X; 

% Potentiel
V = -1 + 2*rand(N,1);

% Fonctions d'interpolation
signal = @(r) bessely(0,k*r) + 1i*besselj(0,k*r);
green  = @(r) signal(r);
d0quad = @(r,rho) besselj(0,r*rho);

% Graphique
% figure
% plot(X(:,1),X(:,2),'*r')
% hold on
% plot(Y(:,1),Y(:,2),'*b')
% hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALCUL DIRECT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% Initialisation
Nspl   = 100;
indBEM = randsample(N,Nspl);
MVref  = zeros(size(indBEM));

% Boucle par ligne
for n = 1:Nspl
    % Distance
    rxy = sqrt( ...
        (X(indBEM(n),1) - Y(:,1)).^2 + ...
        (X(indBEM(n),2) - Y(:,2)).^2 );
    
     % Noyau   
     Gr = green(rxy);
     Gr(rxy<1e-8) = 0 + 1i;
     MVref(n) = Gr.' * V;
end
disp(['Full product      (s) : ',num2str(toc)])
disp(' ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRECOMPUTATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Quadrature 1D 
tic

% Dimensions caracteristiques
dist = @(x,x0) max( sqrt( (x(:,1)-x0(1)).^2 + (x(:,2)-x0(2)).^2) );
rX   = dist(X,mean(X));
rY   = dist(Y,mean(Y)); 
rXY  = dist(mean(X),mean(Y));

% Distance Max
rMax = rXY + rX  + rY;

% Distance SCSD
rMin = rXY - rX - rY;
rMin = max( rMin , 2*rMax*N.^(-1/2) );

% Descente de Newton sur la fonction frequentielle pour les zeros de Y0
f  = @(r) bessely(0,r);
df = @(r) - bessely(1,r);
r1 = (1:1e4)*pi-3*pi/4;
r1 = newton(f,df,r1);
if max(abs(f(r1))) > 1e-12
    BahNonGros
end

% Taille de l'intervalle
L = r1(find(r1>k*rMax,1,'first'))/k;
if abs(real(signal(L))) > 1e-12
    BahNonGros
end

% Descente de Newton sur la fonction frequentielle pour les zeros de J0
f  = @(r) besselj(0,r);
df  = @(r) - besselj(1,r);
r0 = (1:1e4)*pi - pi/2;
r0 = newton(f,df,r0);
if max(abs(f(r0))) > 1e-12
    BahNonGros
end
irho = find(r0>k*rMax,1,'first');

% Discretisation radiale
rQuad = (rMin:(L-rMin)*3e-4:L)';
rSol  = (rMin:(L-rMin)*3e-5:L)';

% Signal avec prolongement impair a gauche, libre a droite
RHS = signal(rQuad);
ref = signal(rSol);

% Calcul du LHS
Nrho  = ceil(sqrt(N));
Mquad = zeros(length(rQuad),Nrho);
Msol  = zeros(length(rSol),Nrho);

% Composante constante
rho(1)     = k;
Mquad(:,1) = d0quad(rQuad,rho);
Msol(:,1)  = d0quad(rSol,rho);

% Construction par frequence
n = 1;
while 1
    % Frequences positives symetriques autour de k
    if irho>n
        rho = [rho , r0(-n+irho)/L , r0(n-1+irho)/L ];
        Nrho = length(rho);  
        Mquad(:,Nrho-1:Nrho) = d0quad(rQuad,rho(end-1:end));
        Msol(:, Nrho-1:Nrho) = d0quad(rSol,rho(end-1:end));
        
    else
        rho = [rho , r0(n-1+irho)/L ];
        Nrho = length(rho);  
        Mquad(:,Nrho) = d0quad(rQuad,rho(end));
        Msol(:, Nrho) = d0quad(rSol,rho(end));
    end

    % Resolution systeme lineaire
    w0 = Mquad(:,1:Nrho) \ RHS;
    
    % Numerical solutions
    sol = Msol(:,1:Nrho) * w0;
    
    % Conditions de sortie de boucle
    if norm(ref-sol,'inf')/norm(ref,'inf') < tol%max(abs(ref-sol))/abs(signal(rMin)) < tol
        break
    end
   
    % Incrementation
    n = n + 1;
end

% Tri par ordre croissant
[rho,ind] = sort(rho);
w0 = w0(ind);

% Infos
disp(['Quadrature 1D     (#) : ',num2str(Nrho)])
disp(['Elapsed time      (s) : ',num2str(toc)])
figure
subplot(2,2,1)
r = (0:1e-3:3*L)'; 
r = r(abs(r)>0.01);
plot(r , signal(r) , 'b' , r, (d0quad(r,rho) * w0) , ' r')
grid on; xlabel('r'); ylabel('signal');
subplot(2,2,2)
plot(rho,w0,'+-r')
grid on; xlabel('rho'); ylabel('weights');
subplot(2,2,3:4)
plot(rSol,(ref-sol)./norm(ref,'inf'))
grid on; xlabel('r'); ylabel('erreur relative L1');


%%% Quadrature 3D en Fourier par spheres
tic

% Initialisation
Np = 0;
Xi  = cell(length(rho),1); wXi = Xi; 

% Construction par frequence
for i = 1:length(rho)
    % Initialisation
    err = 1e6;
    
    while err > tol
        % Incrementation de la discretisation
        Np = Np + 2;
        
        % Discretisation circulaire
        theta = (1:Np)' * (2*pi)/Np ;
        
        % Points d'integration
        S1 = [cos(theta) , sin(theta)];
        
        % Poids d'integration
        wS1 = 2*pi/Np * ones(size(theta));
        
        % Calcul de J0(k*r) = 1/(2*pi) \int_S1 exp(1i*k*r*(s.ej)) pour ej base 2D
        rTest = rho(i)*rMax+pi/2;
        ref = besselj(0,rTest);
        sol = [0 0];
        for j = 1:2
            sol(j) = 1/(2*pi) * (wS1' * cos(rTest * S1(:,j)));
        end
        
        % Erreur
        err = max(abs(ref-sol))/abs(ref);
    end
    
    % Points de quadrature
    Xi{i} = rho(i) * S1;
    
    % Poids de quadrature
    wXi{i} = w0(i) * 1/(2*pi) .* wS1;
end

% Conversion matricielle
Xi  = cell2mat(Xi);
wXi = cell2mat(wXi);
Nxi = length(wXi);

% Infos
disp(['Quadrature 3D     (#) : ',num2str(Nxi)])
disp(['Elapsed time      (s) : ',num2str(toc)])

%%% Close interactions
tic

% Recherche des interactions proches, telles que |x-y| < rMin
[I,rxy] = rangesearch(Y,X,rMin);
jdx = cell2mat(I')';
rxy = cell2mat(rxy')';
idx = zeros(size(jdx));
j = 1;
for i=1:length(I);
    idx(j:j+length(I{i})-1) = i;
    j = j + length(I{i});
end

% Evaluation radiale des termes correctifs
n   = 1e4;
r   = (1.1*rMin) * (1/n:1/n:1)';
Gr = 1 - (d0quad(r,rho) * w0)./signal(r);

% Matrice corrective
Gxy = interp1(r,Gr,rxy,'spline') .* green(rxy);
Gxy(rxy < 1e-8) = 1i - sum(w0);
M = sparse(idx,jdx,Gxy,N,N);

% Infos
disp(['Close corrections (#) : ',num2str(nnz(M))])
disp(['Elapsed time      (s) : ',num2str(toc)])
disp(' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATRIX-VECTOR PRODUCT %%%%%%%%%%%%%%%%%%%%%%%%
tic

% ESPACE Y => FOURRIER Xi
% Vhat = exp(-1i*Xi*Y')*V;
Vhat = nufft2d3(size(Y,1), Y(:,1), Y(:,2), ...
        V, -1, tol, Nxi, Xi(:,1), Xi(:,2) );

% INTEGRATION EN FOURIER
Vhat = wXi .* Vhat;

% FOURIER => X SPACE
% MV = M*V + exp(1i*X*Xi') * Vhat;
MV = M*V + nufft2d3(Nxi, Xi(:,1), Xi(:,2), ...
    Vhat, +1, tol, size(X,1), X(:,1), X(:,2));

% FINAL ERROR
errBEM = norm(MVref-MV(indBEM),'inf')/norm(MVref,'inf');

% CLOSE CORRECTIONS
disp(['Close constant        : ',num2str(nnz(M)/N)])
disp(['Far constant          : ',num2str(Nxi/N)])
disp(['MV product      (INF) : ',num2str(errBEM,'%3.2e')])
disp(['Elapsed time      (s) : ',num2str(toc)])
disp(' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('===> Et voila !')
