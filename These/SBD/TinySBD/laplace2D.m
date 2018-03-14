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
clear all %#ok
close all
clc
alinea = '   ';

% Initialisation
addpath('libGgNufft2D')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETRES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dimension du probleme
N   = 1e5;
tol = 1e-3;

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
signal = @(r) log(r);
green  = @(r) signal(r) ;
d0quad = @(r,rho) besselj(0,r*rho);

% Graphique
% figure
% plot(X(:,1),X(:,2),'*r')
% hold on
% plot(Y(:,1),Y(:,2),'*b')
% hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% CALCUL DIRECT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('I°) VALIDATION ...')

tBEM = tic;

% Initialisation
Nspl   = fix(1*10^8/N)+1;
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
     Gr(rxy<1e-8) = 0;
     MVref(n) = Gr.' * V;
end
clear rxy
clear Gr
disp(['Validation product     (s) : ',num2str(toc(tBEM))])
disp('DONE VALIDATING!')
fprintf('\n\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tTotalTime = tic;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% PRECOMPUTATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('II°) PRECOMPUTATION ')
disp([alinea 'a) Fourier quadrature'])
disp([alinea alinea '-Radial'])
tPrecomputation = tic;
%%% Quadrature 1D 
tQuad1D = tic;

% Dimensions caracteristiques
dist = @(x,x0) max( sqrt( (x(:,1)-x0(1)).^2 + (x(:,2)-x0(2)).^2) );
rX   = dist(X,mean(X));
rY   = dist(Y,mean(Y)); 
rXY  = dist(mean(X),mean(Y));

% Distance Max
rMax = rXY + rX  + rY;

% Distance SCSD
rMin = rXY - rX - rY;
rMin = max( rMin ,0.8*rMax*N.^(-1/2) );

% Discretisation radiale
rQuad = (rMin:(rMax-rMin)*1e-3:rMax)';
rSol  = (rMin:(rMax-rMin)*1e-4:rMax)';

% Signal à approximer
ref = signal(rSol);

% Calcul de l'approximation radiale 
tolH1 = tol*norm(ref,'inf');
tolInf = tol*norm(ref,'inf');
a = rMin/rMax; % Echelle réduite
b = 1; % Dans le cas du log
func = signal; 
derivative = @(x)(1./x); % Dérivée du log
askGraph = false; % Ne pas demander le graphe
% [~,rho,w0] = computeBesselCoeffH1(a,b,func,derivative,tolH1,tolInf,askGraph);
[~,rho,w0] = computeBesselCoeffH1_bis(a,b,func,derivative,tolInf,askGraph);
Nrho = length(rho); % Nombre de fréquences de la quadrature

% Passage à l'échelle
rho = [0 , rho']/rMax;
w0  = [signal(rMax) ; w0];

% Solution à comparer à la référence
% Msol = besselj(0,rSol*rho);
% sol = Msol*w0; 


% Infos
tQuad1D = toc(tQuad1D);
disp([alinea alinea 'Radial quadrature         (#) : ',num2str(Nrho)])
disp([alinea alinea 'Elapsed time              (s) : ',num2str(tQuad1D)])
fprintf('\n')
% Affichage
% figure
% subplot(2,2,1)
% r = (-rMax:1e-3:rMax)'; 
% r = r(abs(r)>0.01);
% plot(r , signal(abs(r)) , 'b')
% hold on
% plot(r, (d0quad(abs(r),rho) * w0) , ' r')
% hold off;
% grid on; xlabel('r'); ylabel('signal');
% subplot(2,2,2)
% plot(rho,w0,'+-r')
% grid on; xlabel('rho'); ylabel('weights');
% subplot(2,2,3:4)
% plot(rSol,(ref-sol)./norm(ref,'inf'))
% grid on; xlabel('r'); ylabel('Relative error Linf');

% clear Msol
% clear sol

%%% Quadrature 2D en Fourier par cercles

disp([alinea alinea '-Circular'])
tQuad2D = tic;

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

tQuad2D = toc(tQuad2D);
disp([alinea alinea 'Circular quadrature       (#) : ',num2str(Nxi)])
disp([alinea alinea 'Elapsed time              (s) : ',num2str(tQuad2D)])
fprintf('\n')
%%% Close interactions
disp([alinea 'b) Close interactions...'])
disp([alinea alinea '-Range search'])
tRangeSearch = tic;

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
% Save memory 
clear I;


tRangeSearch = toc(tRangeSearch);
disp([alinea alinea 'Elapsed time (s)             : ',num2str(tRangeSearch)])
disp([alinea alinea 'Number of close interactions : ',num2str(length(rxy))])
fprintf('\n')

% Evaluation radiale des termes correctifs
disp([alinea alinea '-Assembling'])
tAssembling = tic;
n   = 1e4;
r   = (1.1*rMin) * (1/n:1/n:1)';
Gr = 1 - (d0quad(r,rho) * w0)./signal(r);

% Matrice corrective
Gxy = interp1(r,Gr,rxy,'spline') .* green(rxy);
Gxy(rxy < 1e-8) = 0 - sum(w0);
M = sparse(idx,jdx,Gxy,N,N);
clear rxy; % save memory !!!
clear idx;
clear jdx;

% Infos
tAssembling = toc(tAssembling);
disp([alinea alinea 'Elapsed time      (s)        : ',num2str(tAssembling)])
disp(' ')

disp('DONE PRECOMPUTING ! ')
fprintf('Total precomputing time (s)   : %s \n\n',num2str(toc(tPrecomputation)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATRIX-VECTOR PRODUCT %%%%%%%%%%%%%%%%%%%%%%%%
disp('III°) MATRIX VECTOR PRODUCT :')
disp([alinea 'a) Far field (NuFFT)'])
tMV = tic;
tNUFFT = tic;

% ESPACE Y => FOURRIER Xi
% Vhat = exp(-1i*Xi*Y')*V;
Vhat = nufft2d3(size(Y,1), Y(:,1), Y(:,2), ...
        V, -1, tol, Nxi, Xi(:,1), Xi(:,2) );

% INTEGRATION EN FOURIER
Vhat = wXi .* Vhat;

% FOURIER => X SPACE
% MV = M*V + exp(1i*X*Xi') * Vhat;
MV = nufft2d3(Nxi, Xi(:,1), Xi(:,2), ...
    Vhat, +1, tol, size(X,1), X(:,1), X(:,2));

% Infos
tNUFFT = toc(tNUFFT);
disp([alinea 'Elapsed time      (s) : ',num2str(tNUFFT)])
disp(' ')

% Close corrections : 
disp([alinea  'b) Close correction (Sparse product)'])
tCloseCorrec = tic;
MV = MV + M*V;

tCloseCorrec = toc(tCloseCorrec);
disp([alinea 'Elapsed time      (s) : ',num2str(tCloseCorrec)])
disp(' ')

disp('DONE COMPUTING MATRIX-VECTOR PRODUCT ! ')
tMV = toc(tMV);
fprintf('Total matrix-vector time(s)   : %s \n\n',num2str(tMV));


% FINAL ERROR
errBEM = norm(MVref-MV(indBEM),'inf')/norm(MVref,'inf');

% SUMMARY 
fprintf('\n\n')
disp('SUMMARY')
disp(['Close constant        : ',num2str(nnz(M)/N)])    
disp(['Far constant          : ',num2str(Nxi/N)])
disp(['Error                 : ',num2str(errBEM,'%3.2e')])
disp(['Total time (s)        : ',num2str(toc(tTotalTime))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Graphic for computing time 
figure
pieX = ([tQuad1D, tQuad2D, tRangeSearch, tAssembling, tNUFFT, tCloseCorrec]);
pieLabs = {'Quad 1D', 'Quad 2D', 'Range search', 'Assembling','NUFFT','MV sparse'};

pie(pieX,pieLabs);
