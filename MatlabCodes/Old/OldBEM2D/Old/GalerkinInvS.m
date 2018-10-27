function [ lambda,K ] = GalerkinInvS(f,X,Pk)
% L'objet de cette fonction est d'obtenir la solution de la formulation
% variationnelle discrète de l'équation S lambda = f où S est le potentiel
% de simple couche, sur le segment (-1,1). 
% Le pas du maillage est h = pi/(N-1). N est le nombre de noeuds du
% maillage, il y a N-1 mailles, et N-1 points d'intégrations. 
% Le vecteur lambda obtenu représente les valeurs prises par la fonction
% contante par morceaux lambda sur chaque maille. 

K = rigidite_simpleCouche(X,Pk);
L = secMembre(f,X);
L = L(:);
lambda = K\L;



end

