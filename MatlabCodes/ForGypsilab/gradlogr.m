function [compNormale,compTangentielle] = gradlogr(A,B,X,N,tau)
% Cette fonction calcule, pour un point X dans le plan, et un segment
% [A,B], la valeur principale de l'intégrale
% \int_{[A,B]} \frac{(Y - X) dY}{\norm{Y - X}^2}
% C'est-à-dire p.v \int_{[A,B]} \grad_y G_0(Y-X) dY
% Inputs : 
% A,B   : les coordonnées des extrémités du segement
% N     : vecteur normal unitaire au segment [A,B] orienté d'un côté. 
% X     : vecteur de points auxquels on souhaite calculer l'intégrale.
%       La composante normale devient opposée si on change N en -N
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
compNormale = -atan(b./d) + atan(a./d);
% Pour d = 0, l'intégrale vaut 0. (C'est faux si d est non nul, même très
% proche de 0, il y a une discontinuité). 
compNormale(d == 0) = 0;

% Composante tangentielle (valeur principale).
% Si d est non nul, l'intégrale est convergente :
compTangentielle = 1/2*log((b.^2 + d.^2)./(a.^2 + d.^2));
% Si d est nul, on prend la valeur principale
compTangentielle(d==0) = 1/2*log(abs(b(d==0)./a(d==0)));

    function[out] = scal_vec(X,n)
        % Calcule le tableau T tq T(i) = X_i scalaire n
        % X tableau N x 2 : N points dans R^2
        % n vecteur 1 x 2 : un vecteur de R^2
        
        out = X(:,1)*n(1) + X(:,2)*n(2);
        
    end





end

