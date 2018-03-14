function [ K ] = rigidite_simpleCoucheP0( X )

% Les points du maillage sont donnés dans la variable X, de taille
N = length(X);
% On calcule les points du maillage après changement de variable 
T = acos(X);
% Du fait que Arccos est décroissante, on remet les nouveaux points d'intégration à
% l'endroit. 
T = T(end:-1:1);

% Le calcul de la matrice de rigidité est décomposé en 4 parties : A, B, C et
% D.
A = matA(T,N);
B = matB(T,N);
C = rot90(rot90(B,1),1); % fliplr(flipud(B));
D = matD(T,N);

K = -(A + B + C + D)/(2*pi);

end

function[A] = matA(T,N)


A = zeros(N-1,N-1);

f = @(c)(c^2*log(c));
g = @(c,d)(-3/2*c*d);
h = @(i,j)(f(T(j) - T(i)));
l = @(i,j)(g((T(i+1)-T(i)),(T(j+1)-T(j))));

N = length(T);
for i = 1:N-1
    for j = i:N-1
        if j==i
            A(i,i) = h(i,i+1) + l(i,i);
        elseif j==i+1
            A(i,i+1) = 1/2*h(i,i+2) - 1/2*h(i,i+1) - 1/2*h(i+1,i+2) + l(i,i+1);
        else
            A(i,j) = 1/2*h(i,j+1) - 1/2*h(i,j) - 1/2*h(i+1,j+1) + 1/2*h(i+1,j) + l(i,j);
        end
    end
end

A = symmetrize(A);

end


function[B] = matB(T,N)

B = zeros(N-1,N-1);

f = @(c)(c^2*log(c));
g = @(c,d)(-(3/2+log(2))*c*d);
h = @(i,j)(g((T(i+1)-T(i)),(T(j+1)-T(j))));

for i = 1:N-1
    for j = i:N-1
        if and(T(i)==0,T(j)==0);
            B(i,i) = 1/2*(T(i+1)-T(i))^2*(-3 + 2*log(T(i+1) - T(i))+2*log(2));
        else
            B(i,j) =    1/2*f(  T(i+1) +  T(j+1))  -  1/2*f(T(i)+ T(j+1)) - ...
                        1/2*f(  T(i+1)   +  T(j))  +  1/2*f(T(i) + T(j))...
                        + h(i,j);
        end
    end
end

B = symmetrize(B);

end

function[D] = matD(T,N)
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

% points d'intégration
Delta = diff(T);
phi = T(1:end-1) + Delta/2; 

phi_i = repmat(reshape(phi,1,N-1),N-1,1);
phi_j = repmat(reshape(phi,N-1,1),1,N-1);

Delta_i = repmat(reshape(Delta,1,N-1),N-1,1);
Delta_j = repmat(reshape(Delta,N-1,1),1,N-1);

f1 = @(u)(sinc(u/pi)); % attention à la déf de sinc sous matlab...
f2 = @(u)((u<=pi/2).*sinc(u/pi)./(pi-u) + (u>pi/2).*sinc((pi-u)/pi)./u);
integrande = @(t1,t2)(log(f2((t1+t2)/2).*f1(abs(t1-t2)/2)));

D = Delta_i.*Delta_j.*integrande(phi_i,phi_j);
end




