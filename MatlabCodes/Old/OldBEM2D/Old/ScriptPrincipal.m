%% Résolution du problème de Dirichlet sur le segment (-1,1)
clear all; %#ok
close all;
clc;

Main;

% Nombre de points du maillage
N = 50;


% Donnée de Dirichlet
n = 7;
if n==0
    lambda_n = log(2)/2;
else
    lambda_n = 1/(2*n);
end
f = @(t)(chebyshevT(n,t)*lambda_n);

% La variable d'espace sera notée x. Le nombre k est le nombre d'onde. 

% On souhaite résoudre le problème Laplace avec conditions de Dirichlet dans R²\(-1,1) avec u =
% uinc sur (-1,1) et condition de radiation de Sommerfeld. On se ramène à
% une équation intégrale sur la frontière du domaine : 
%
%           S[l] = uinc 
% Avec S le simple couche, [l] le saut de dérivée normale. Le champ rayonné est ensuite obtenu
% via S[l] partout dans l'espace. 
% On maille la frontière du domaine de manière à prendre en compte la
% singularité, en faisant le changement de variable circulaire [l] = a/w
% où 

w = @(t)(sqrt(1-t.^2)); 

% t représente un point courant sur (-1,1), càd un x = (t,0) avec -1 < t <
% 1
% On fait un maillage du segment (-1,1) 
plot([-1 1],[0 0]);
axis equal
hold on; 
% formé de N noeuds de Thcebichev, et donc de 
Nmaille = N-1;
% mailles. Les noeuds de Tchebitchev sont donnés par
tNodes = tchebNodes(N);
% il y a 
assert(length(tNodes)==N);
% Noeuds dans ce maillage, qui sont les projetés sur (-1,1) d'un maillage
% uniforme du cercle unité, comme le montre ce graphe

cercleUnite = [cos((0:500)'*2*pi/500), sin((0:500)'*2*pi/500)];
plot(cercleUnite(:,1),cercleUnite(:,2)); 
plot(tNodes',w(tNodes'),'r *');
plot(tNodes',-w(tNodes'),'g *');
plot(tNodes,0*tNodes,'k *');
for j = 1:N
    plot([tNodes(j), tNodes(j)],[0,w(tNodes(j))],'k--');
end
title('Maillage de Tchebitchev')


% On résout le problème de manière variationnelle, et on 
% se ramène à l'inversion d'un système linéaire donné par 
%
%                   K U = L
%
% Où 
K = matriceDeRigidite(Nmaille);
% On peut voir à quoi ressemble cette matrice : 
hold off
figure
axis xy
imagesc(K);
title('Matrice de rigidité')

% Et où 
L = secondMembre(f,Nmaille);

% Ainsi, 
U = K\L; 
% Donne une fonction constante par morceaux, valant U(i) sur chaque maille,
% qui représente le saut de la dérivée normale. On va pouvoir utiliser la
% toolbox 2D que j'avais faite pour calculer le champ rayonné. 

figure
plot(tNodes(1:end-1) + diff(tNodes)/2,U,'x');
hold on
plot(tNodes,chebyshevT(n,tNodes));
title('Exact vs. Approx')















