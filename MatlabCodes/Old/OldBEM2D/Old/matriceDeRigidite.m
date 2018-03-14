function[K] = matriceDeRigidite(Nmailles)

% La matrice exacte à calculer est
%
%       Kij = \int_{x_i}^{x_{i+1}} \int_{x_j}^{x_{j+1}} \frac{\ln(|x - x'|){\sqrt{(1-|x|^2)(1-|x'|^2)}}dxdx'
%
% Où les x_j sont les neouds de Tchebitchev
% Après changement de variable, on obtient
%
%       Kij = \int_{\theta_i}^{\theta_{i+1}} \int_{\theta_j}^{\theta_{j+1}} \frac{\ln(|cos(\theta) - \cos(\theta')|)d\theta d\theta'
%
% Où l'integrande est devenue intégrable.
N = Nmailles + 1; % Le nombre de points de Tchebitchev vaut le nombre de mailles + 1
%
% Pour obtenir une approximation de l'intégrale, on la découpe en 4
% morceaux, selon le document
% /home/martin/Documents/These/ComptesRendusDeThese/RigiditeSimpleCoucheSegment.pdf

% Soit
Delta = pi/(Nmailles);
% L'écart entre deux angles theta_i consécutifs. 

% On remplit désromais la matrice K avec des valeurs approchées :
A = computePartA(Nmailles,Delta);
B = computePartB(Nmailles,Delta);
C = rot90(rot90(B,1),1); % fliplr(flipud(B));
D = computePartD(Nmailles,Delta);
K = -(A+B+C+D)/(2*pi);
end


function[A] = computePartA(Nmailles,Delta)
% Cette partie correspond au calcul de l'intégrale 
%
%       A_{ij} =
%       \int_{\theta_i}^{\theta_{i+1}}\int_{\theta_j}^{\theta_{j+1}}
%       \ln|\theta - \theta'|d\theta d\theta'  
% 
%

I=repmat(0:Nmailles-1,Nmailles,1);
J = I';
A = Delta^2*(log(Delta)+(log(abs(I-J))));
for i=0:(Nmailles-1)
    A(i+1,i+1) = Delta^2*(log(Delta)-3/2);
    if i>0
        A(i,i+1) = Delta^2*(log(Delta) + 2*log(2) - 3/2);
    end
    if i<Nmailles-1
        A(i+2,i+1) = Delta^2*(log(Delta) + 2*log(2) - 3/2);        
    end
end
% On doit avoir
assert(max(max((abs(A-A')))) < 1e-10);
% Car A est symétrique. 
% Aucune valeur dans A ne doit être infinie
assert(max(max(abs(A))) < Inf);

end

function[B] = computePartB(Nmailles,Delta)
% Cette partie correspond au calcul de l'intégrale 
%
%       B_{ij} =
%       \int_{\theta_i}^{\theta_{i+1}}\int_{\theta_j}^{\theta_{j+1}}
%       \ln\frac{\theta + \theta'}{2}d\theta d\theta'  
% 
%
I=repmat(0:Nmailles-1,Nmailles,1);
J=I';
B = Delta^2*(log(Delta)+(log(abs(I+J+1)/2)));
B(1,1) = Delta^2*(log(Delta)+log(2) - 3/2);
% On doit avoir là encore
assert(max((max(abs(B-B')))) < 1e-10);
% Car B est symétrique. 
% Aucune valeur dans B ne doit être infinie
assert(max(max(abs(B))) < Inf)

end

function[D] = computePartD(Nmailles,Delta)
% Cette partie correspond au calcul de l'intégrale 
%
%       D_{ij} =
%       \int_{\theta_i}^{\theta_{i+1}}\int_{\theta_j}^{\theta_{j+1}} 
%       \ln\left(\frac{\sin\frac{\theta + \theta'}{2}}{\frac{\theta +
%       \theta'}{2} \left(\pi-\frac{\theta + \theta'}{2}\right)}\sinc
%       \frac{\theta-\theta'}{2}\right)d\theta d\theta' 
%


% L'intégrande est très régulier, on peut faire une quadrature à un point.
% Les points d'intégrations sont, en posant 
N = Nmailles + 1;
% (Nmailles + 1 angles)
theta = anglesTcheb(N);
% Les angles de Tchebitchev
phi = theta(1:end-1) + diff(theta)/2;
% Les milieux des segments theta_i,theta_i+1, et on a bien 
assert(length(phi)==Nmailles);
% Donc un point d'intégration par segment. 

phi_i = repmat(reshape(phi,1,Nmailles),Nmailles,1);
phi_j = repmat(reshape(phi,Nmailles,1),1,Nmailles);

f1 = @(u)(sinc(u/pi)); % attention à la déf de sinc sous matlab...
f2 = @(u)((u<=pi/2).*sinc(u/pi)./(pi-u) + (u>pi/2).*sinc((pi-u)/pi)./u);
integrande = @(t1,t2)(log(f2((t1+t2)/2).*f1(abs(t1-t2)/2)));

D = Delta^2*integrande(phi_i,phi_j);

assert(max(max(abs(D-D')))<1e-10);
end





