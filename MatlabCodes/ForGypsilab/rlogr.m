function [compNormale,compTangentielle] = rlogr(A,B,X,N,tau)
% Cette fonction calcule, pour un point X dans le plan, et un segment
% [A,B], la valeur de l'intégrale
% \int_{[A,B]} (Y - X) log( || Y - X ||) dY
% Inputs :
% A,B   : les coordonnées des extrémités du segement
% X     : vecteur de points auxquels on souhaite calculer l'intégrale.
% N     : vecteur normal unitaire au segment [A,B] orienté d'un côté.
%       La composante normale de l'intégrale devient bien sûr opposée si on
%       change N en -N
% tau   : vecteur tangentiel unitaire du segment [A,B], orienté de A vers
%       B.

% Paramétrisation de l'intégrale

XA = A - X;
a = scal_vec(XA,tau);
d = -scal_vec(XA,N);
XB = B - X;
b = scal_vec(XB,tau);

% Composante normale.
% Pour d non nul, l'intégrale vaut arctan(b/d) - arctan(a/d)
F = @(x)(1/2*x.*log(d.^2 + x.^2)-x + d.*atan(x./d));
compNormale = d.*(F(b) - F(a));
% Pour d = 0, l'intégrale vaut 0. (C'est faux si d est non nul, même très
% proche de 0, il y a une discontinuité).
F = @(x)(xlog(x) - x);
compNormale(d == 0) = F(b(d == 0)) - F(a(d == 0));

% Composante tangentielle (valeur principale).
% Si d est non nul, l'intégrale est convergente :
F = @(x)(1/4*(xlog(d.^2 + x.^2) - x.^2));
compTangentielle = F(b) - F(a);
    function[out] = scal_vec(X,n)
        % Calcule le tableau T tq T(i) = X_i scalaire n
        % X tableau N x 2 : N points dans R^2
        % n vecteur 1 x 2 : un vecteur de R^2
        
        out = X(:,1)*n(1) + X(:,2)*n(2);
        
    end
    function[out] = xlog(x)
        % Calcule le tableau T tq T(i) = X_i scalaire n
        % X tableau N x 2 : N points dans R^2
        % n vecteur 1 x 2 : un vecteur de R^2
        
        out = x.*log(x);
        out(abs(x) < 1e-12) = 0;
        
    end





end

