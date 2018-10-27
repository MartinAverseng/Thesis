function[L] = secondMembre(f,Nmaille)
% Le vecteur à calculer
% 
%           I_i = \int_{x_i}^{x_{i+1}} f(x)/\sqrt{1-x^2}dx
% 
% Après changement de variable, on obtient 
% 
%           I_i = \int_{\theta_i}^{\theta_{i+1}} f(\cos\theta)d\theta 
% 
% On utilise alors une formule à un point 
% On commence par définir les paramètres suivants, voir matriceDeRigidite.m
N = Nmaille + 1;
theta = anglesTcheb(N);
phi = theta(1:end-1) + diff(theta)/2;
Delta = pi/(N-1);
% Le vecteur L est alors tout simplement

L = Delta*f(cos(phi));

end