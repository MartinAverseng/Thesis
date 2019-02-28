% Test de l'idée de l'initégrale à poids sur un triangle ......

clear all
close all

rhos = rand(1,3);
phis = 2*pi*(rand(1,3)-1/2);

ABC = [rhos.*cos(phis); rhos.*sin(phis)];

A = ABC(:,1);
B = ABC(:,2);
C = ABC(:,3);

% Plot the triangle

plot([A(1),B(1)],[A(2),B(2)]);
hold on
plot([B(1),C(1)],[B(2),C(2)]);
plot([A(1),C(1)],[A(2),C(2)]);

% Détermination des variableflecheC = AB/2 + nC;s phi_A, phi_B, phi_C

phis_guess = atan2(ABC(2,:),ABC(1,:));

% Normales orientées aléatoirement pour commencer :
rotation = [0 -1;1,0];
AB = B-A; nC = rotation*AB;
BC = C-B; nA = rotation*BC;
CA = A-C; nB = rotation*CA;

if sum(AB.*nA) < 0
    nA = -nA;
end

if sum(BC.*nB) < 0
    nB = -nB;
end

if sum(CA.*nC) < 0
    nC = -nC;
end

quiver(A(1)/2 + B(1)/2,A(2)/2 + B(2)/2,nC(1),nC(2));
quiver(B(1)/2 + C(1)/2,B(2)/2 + C(2)/2,nA(1),nA(2));
quiver(C(1)/2 + A(1)/2,C(2)/2 + A(2)/2,nB(1),nB(2));


% Test if O is in the triangle: 
inside = and(and(sum(B.*nC)>=0,sum(C.*nA)>=0),sum(A.*nB)>=0);
disp(inside);

plot(0,0,'*');
axis equal


% If it is inside the triangle, we integrate in each of the three inner
% triangles OAB, OBC and OCA. 

% Triangle OAB

phiA = phis_guess(1);
phiB = phis_guess(2);

[phim,argm] = min([phiA,phiB]);
[phiM,argM] = max([phiA,phiB]);

rhoAB = [rhos(1),rhos(2)];
rho0 = rhoAB(argm);

% Gauss and quadrature points in phi
NquadPhi = 4;
Nquadrho = 3;
[phi,wphi] = Gauss_Legendre1D(NquadPhi,phim,phiM);
for i = 1:NquadPhi
    % Determine rhoMax(phi)
    rhoMax = sum(A.*nC)/sum([cos(phi(i));sin(phi(i))].*nC);
    plot([0,rhoMax*cos(phi(i))],[0,rhoMax*sin(phi(i))],'--');
    % Build the quadrature in rho: 
    [theta_i,wi] = Gauss_Legendre1D(Nquadrho,acos(rhoMax),pi/2);
    wphi_rho = wphi(i)*(wi.*cos(theta_i));
    rho_i = cos(theta_i);
    plot(rho_i*cos(phi(i)),rho_i*sin(phi(i)),'x')
end




