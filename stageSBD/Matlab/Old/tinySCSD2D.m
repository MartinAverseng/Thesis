function tinySCSD2D()
% function tinySCSD()
% Authors: matthieu.aussal@polytechnique.edu
% Creation : 16.12.2015
% Copyright © Ecole polytechnique, Matthieu AUSSAL, 2015.

% Nettoyage ecran
clear all
close all
clc

% Nufft from L. Greengard
addpath('libGgNufft2D')

% Dimension du probleme
N    = 1e3;
tol  = 1e-3;

% Nuages de points dansl'espace 3D
X = -1 + 2*rand(N,3);
Y = -1 + 2*rand(N,3);

% Dimensions caracteristiques des interactions
rMax  = 2*sqrt(3);
rMin  = rMax/10;

% Potentiel
V = -1 + 2*rand(N,1);

% Noyau de green
green1D = @(r)(-1/(2*pi)*log(r));
green3D = @(r,x,y) (-1/(2*pi)*log(r));
quad = @(r,rho) besselj(0,r*rho);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FOR LOOPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% Boucle for
MVref = zeros(N,1);
for i = 1:N
    rxy = sqrt( ...
        (X(i,1) - Y(:,1)).^2 + ...
        (X(i,2) - Y(:,2)).^2 + ...
        (X(i,3) - Y(:,3)).^2 );
    Gkr = green3D(rxy,X(i,:),Y);
    Gkr(rxy<1e-8) = 0;
    MVref(i) = Gkr.' * V;
end
disp(['Full product           (s) : ',num2str(toc)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRECOMPUTATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% Discretisation of radius for mono dimensional kernel (fct of r)
rSol  = (rMin:1e-3*rMin:rMax)';

% Approximation en bessel serie
func = @(r)(green1D(r*(rMax+rMin)));
weights = besselQuad(func,rMin/(rMax+rMin),tol);
Nrho = length(weights);
rho = besselZeros(Nrho)'/(rMax+rMin);

% Numerical solutions
sol = quad(rSol,rho) * weights;
ref = green1D(rSol);


disp(['SCSD 1D quadrature     (s) : ',num2str(toc)])

% Graphique
figure(1)
plot(rSol,ref,'b')
hold on
plot(rSol,sol,'r--')
figure(2)
plot(rho,weights,'+-')

% All close interactions for |x-y| < rMin
tic
[I,rxy] = rangesearch(Y,X,rMin);
jdx = cell2mat(I')';
rxy = cell2mat(rxy')';
idx = zeros(size(jdx));
j = 1;
for i=1:length(I);
    idx(j:j+length(I{i})-1) = i;
    j = j + length(I{i});
end

% Exact value of the green kernel
Gxy = green1D(rxy);
Gxy(rxy < 1e-8) = 0;
% SCSD values
SCSDxy = quad(rxy,rho) * weights;

% Matrix value
val = Gxy - SCSDxy;

% Compute correction matrix
M = sparse(idx,jdx,val);
disp(['SCSD corrective matrix (s) : ',num2str(toc)])

% Graphique
figure(3)
spy(M)

% 2D fourier quadrature
tic
Xi  = cell(length(rho),1); wXi = Xi; Nang = 1;
for j = 1:length(rho)
    % Quadrature spherique de rayon k*rho, discretisee par Nang
    [Srho,wSrho,Nang] = oprSphericalQuadrature(Nang,rho(j)*rMax,tol);
    
    % Points de quadrature
    Xi{j} = rho(j) * Srho;
    
    % Poids de quadrature
    wXi{j} = weights(j) * rho(j)/(4*pi) .* wSrho;
end
Xi = cell2mat(Xi);
wXi = cell2mat(wXi);
Nxi = length(wXi);
disp(['SCSD 3D quadrature (s)     : ',num2str(toc)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATRIX-VECTOR PRODUCT %%%%%%%%%%%%%%%%%%%%%%%%
tic

% Y SPACE => FOURRIER
% Vhat = exp(-1i*Xi*Y')*V;
Vhat = nufft3d3(size(Y,1), Y(:,1), Y(:,2), Y(:,3), ...
    V, -1, tol, Nxi, Xi(:,1), Xi(:,2), Xi(:,3));

% FOURIER INTEGRATION
Vhat = wXi .* Vhat;

% FOURIER => X SPACE
% MV = exp(1i*X*Xi') * Vhat;
MV = nufft3d3(Nxi, Xi(:,1), Xi(:,2), Xi(:,3), ...
    Vhat, +1, tol, size(X,1), X(:,1), X(:,2), X(:,3));

% CLOSE CORRECTIONS
MV = MV + M*V;
disp(['SCSD matrix vecor product (s) : ',num2str(toc)])

% Final constants
farCte = Nxi/N
closeCte = nnz(M)/N

% FINAL ERROR
relativeErrorL2 = norm(MVref-MV)/norm(MVref)


disp('===> Et voila !')
end

function [S2,wS2,Nang] = oprSphericalQuadrature(Nang,krMax,tol)
% function tinySCSD()
% Authors: matthieu.aussal@polytechnique.edu
% Creation : 16.12.2015
% Copyright © Ecole polytechnique, Matthieu AUSSAL, 2015.

% Donnees pour le calcul de l'erreur
ref = sin(krMax);
sol = zeros(1,3);

% Verification que |ref| > tol;
while abs(ref) < tol
    krMax = krMax + 0.1;
    ref = sin(krMax);
end

% Constante multiplicative
cte = (krMax)/(4*pi);

% Boucle sur la discretisation angulaire
while 1
    % Sphere gauss Legendre
    [~,phi,theta,wg] = sphereQuad(1,Nang,2*Nang,1);
    phi = phi - pi/2;
    [x,y,z] = sph2cart(theta,phi,1);
    nx = numel(x);
    S2 = [reshape(x,nx,1),reshape(y,nx,1),reshape(z,nx,1)];
    wS2 = 3 * wg;
    
    % Calcul de sin(k*r) = (k*r)/(4*pi) \int_S2 exp(1i*k*r*(s.ej)) pour ej base 3D
    for j = 1:3
        sol(j) = cte * (wS2' * cos(krMax * S2(:,j)));
    end
    
    % Erreur
    err = max(abs(ref - sol))/abs(ref);
    
    % Incrementation
    if err >= tol
        Nang = Nang+1;
    else
        break
    end
end
end


function [r,t,p,w] = sphereQuad(nr,nt,np,rad)
% Written by: Greg von Winckel - 04/13/2006
% Contact: gregvw(at)math(dot)unm(dot)edu
% URL: http://www.math.unm.edu/~gregvw
[r,wr]=rquad(nr,2);         % radial weights and nodes (mapped Jacobi)
if rad==inf                 % infinite radius sphere
    
    wr=wr./(1-r).^4;        % singular map of sphere radius
    r=r./(1-r);
else                        % finite radius sphere
    wr=wr*rad^3;            % Scale sphere radius
    r=r*rad;
end
[x,wt]=rquad(nt,0);
t=acos(2*x-1); wt=2*wt;     % theta weights and nodes (mapped Legendre)
p=2*pi*(0:np-1)'/np;        % phi nodes (Gauss-Fourier)
wp=2*pi*ones(np,1)/np;      % phi weights
[rr,tt,pp]=meshgrid(r,t,p); % Compute the product grid
r=rr(:); t=tt(:); p=pp(:);
w=reshape(reshape(wt*wr',nr*nt,1)*wp',nr*nt*np,1);
end


function [x,w]=rquad(N,k)
% Written by: Greg von Winckel - 04/13/2006
% Contact: gregvw(at)math(dot)unm(dot)edu
% URL: http://www.math.unm.edu/~gregvw
k1=k+1; k2=k+2; n=1:N;  nnk=2*n+k;
A=[k/k2 repmat(k^2,1,N)./(nnk.*(nnk+2))];
n=2:N; nnk=nnk(n);
B1=4*k1/(k2*k2*(k+3)); nk=n+k; nnk2=nnk.*nnk;
B=4*(n.*nk).^2./(nnk2.*nnk2-nnk2);
ab=[A' [(2^k1)/k1; B1; B']]; s=sqrt(ab(2:N,2));
[V,X]=eig(diag(ab(1:N,1),0)+diag(s,-1)+diag(s,1));
[X,I]=sort(diag(X));
x=(X+1)/2; w=(1/2)^(k1)*ab(1,2)*V(1,I)'.^2;
end

