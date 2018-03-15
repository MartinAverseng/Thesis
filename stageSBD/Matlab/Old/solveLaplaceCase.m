function [ tPrecomputation, tMV,   tTotalTime,MV ] = solveLaplaceCase(N,tol,const)


a = const*sqrt(abs(log(tol))/N);
disp(['a = ' num2str(a)]);
alinea = '   ';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETRES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tTotalTime = tic;

% Nuages de points dans un cube + potentiel
X = rand(N,2);
Y = rand(N,2);
V = randn(N,1);

% Noyau de green
signal = @(r) log(r);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%% PRECOMPUTATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('II°) PRECOMPUTATION ')
disp([alinea 'a) Fourier quadrature'])
disp([alinea alinea '-Radial'])
tPrecomputation = tic;
%%% Quadrature 1D

% Dimensions caracteristiques :
% Borne sup pour rMax
dist = @(x,x0) max( sqrt( (x(:,1)-x0(1)).^2 + (x(:,2)-x0(2)).^2) );
rX   = dist(X,mean(X));
rY   = dist(Y,mean(Y));
rXY  = dist(mean(X),mean(Y));
rMax = rXY + rX  + rY;
% rMin
rMin = a*rMax;

% Calcul de l'approximation radiale
tolInf = tol*log(rMax);
a = rMin/rMax; % Echelle réduite
b = 1; % Dans le cas du log
func = signal;
derivative = @(x)(1./x); % Dérivée du log
askGraph = false; % Ne pas demander le graphe
[~,rho,w0] = computeBesselCoeffH1_bis(a,b,func,derivative,tolInf,askGraph);
% Repassage à l'échelle
rho = [0 , rho']/rMax;
w0  = [signal(rMax) ; w0];


%%% Quadrature 2D en Fourier par cercles

disp([alinea alinea '-Circular'])

% Initialisation
Xi  = cell(length(rho),1); wXi = Xi;

% Construction par frequence
gamma3 = 5;
gamma4 = 0.33;
for i = 1:length(rho)
    % Initialisation
    Np = fix(rho(i)*rMax + gamma3 * (rho(i)*rMax)^gamma4)+1;
    % Discretisation circulaire
    theta = (1:Np)' * (2*pi)/Np ;
    % Points d'integration
    S1 = [cos(theta) , sin(theta)];
    % Poids d'integration
    wS1 = 2*pi/Np * ones(size(theta));
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

disp([alinea alinea 'Circular quadrature       (#) : ',num2str(Nxi)])
fprintf('\n')
%%% Close interactions
disp([alinea 'b) Close interactions...'])
disp([alinea alinea '-Range search'])

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
NCI = length(rxy);
disp([alinea alinea 'Number of close interactions : ',num2str(NCI)])
fprintf('\n')

disp([alinea alinea '-Assembling'])

% champ proche (M)
B1 = signal(rxy);

% correction de l'erreur due à la NUFFT (B)
B2 = nufft2d3(Nxi, Xi(:,1), Xi(:,2), ...
    wXi, +1, tol, length(idx), X(idx,1) - Y(jdx,1), X(idx,2) - Y(jdx,2));

% Matrice creuse
B = B1 - B2;
M = sparse(idx,jdx,B,N,N);

clear rxy; % save memory !!!
clear idx;
clear jdx;

% Infos
disp(' ')

tPrecomputation = toc(tPrecomputation);
disp('DONE PRECOMPUTING ! ')
fprintf('Total precomputing time (s)   : %s \n\n',num2str(tPrecomputation));


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

balance = NCI/Nxi;
disp(['Balance close / far : ' num2str(balance)]);
tTotalTime = toc(tTotalTime);

end

