function [X,W] = quadSphTri(A,B,C,N)
% Computes a quadrature on a spherical triangle A B C
% We use the gnomonic projection on the plane tangent to the centroid of
% the triangle. 

if size(A,2)==2
    A(:,3) = real(sqrt(1 - A(:,1).^2 - A(:,2).^2));
    B(:,3) = real(sqrt(1 - B(:,1).^2 - B(:,2).^2));
    C(:,3) = real(sqrt(1 - C(:,1).^2 - C(:,2).^2));   
elseif and(and(A(:,3) == 0, B(:,3) == 0),C(:,3) ==0)
    A(:,3) = real(sqrt(1 - A(:,1).^2 - A(:,2).^2));
    B(:,3) = real(sqrt(1 - B(:,1).^2 - B(:,2).^2));
    C(:,3) = real(sqrt(1 - C(:,1).^2 - C(:,2).^2));   
elseif abs(A(:,1).^2 + A(:,2).^2 + A(:,3).^2 - 1 > 1e-6)
    normA = sqrt(A(:,1).^2 + A(:,2).^2 + A(:,3).^2);
    A(:,1) = A(:,1)./normA; A(:,2) = A(:,2)./normA; A(:,3) = A(:,3)./normA; 
    normB = sqrt(B(:,1).^2 + B(:,2).^2 + B(:,3).^2);
    B(:,1) = B(:,1)./normB; B(:,2) = B(:,2)./normB; B(:,3) = B(:,3)./normB; 
    normC = sqrt(C(:,1).^2 + C(:,2).^2 + C(:,3).^2);
    C(:,1) = C(:,1)./normC; C(:,2) = C(:,2)./normC; C(:,3) = C(:,3)./normC; 
end

O = 1/3*(A + B + C);
normO = sqrt(O(:,1).^2 + O(:,2).^2 + O(:,3).^2);
O(:,1) = O(:,1)./normO; O(:,2) = O(:,2)./normO; O(:,3) = O(:,3)./normO; 

[theta_O,phi_O] = angles(O,0,0);

% The tangent plane is parametrized by the coordinates (t,s) where O =
% (0,0)
% Coordinates of the projections of A,B, and C on this plane:

% sphere
% hold on
% plot3(A(1),A(2),A(3),'.','Markersize',50);
% plot3(B(1),B(2),B(3),'.','Markersize',50);
% plot3(C(1),C(2),C(3),'.','Markersize',50);
% plotOnSphere(theta_O,phi_O);

[thetaA,phiA] = angles(A,theta_O,phi_O);
sA = tan(thetaA).*cos(phiA); tA = tan(thetaA).*sin(phiA);
[thetaB,phiB] = angles(B,theta_O,phi_O);
sB = tan(thetaB).*cos(phiB); tB = tan(thetaB).*sin(phiB);
[thetaC,phiC] = angles(C,theta_O,phi_O);
sC = tan(thetaC).*cos(phiC); tC = tan(thetaC).*sin(phiC);


% Debug : plot the gnomonic projection of the triangle on the tangent plane
% G = @(Y)(Y + (1./dot(O,Y) - 1)*Y);
% Ap = G(A); Bp = G(B); Cp = G(C);

% hold on
% plot3(Ap(1),Ap(2),Ap(3),'.','Markersize',50);
% plot3(Bp(1),Bp(2),Bp(3),'.','Markersize',50);
% plot3(Cp(1),Cp(2),Cp(3),'.','Markersize',50);
% 

% Debug : plot the tangent plane
% [x y] = meshgrid(0:0.1:1); 
% z = -1/O(3)*(O(1)*x + O(2)*y - 1); 
% surf(x,y,z) 


J = @(s,t)(1./(1 + s.^2 + t.^2).^(3/2)); % Jacobian on the plane
theta_x = @(s,t)(atan(sqrt(s.^2 + t.^2)));
phi_x = @(s,t)(atan2(t,s));


% Standard gauss quadrature on the triangle A' B' C'
Ap = [sA tA]; Bp = [sB, tB]; Cp = [sC tC];
xw = TriGaussABC(Ap,Bp,Cp,N);


% Debug
% plot([Ap(1) Bp(1) Cp(1) Ap(1)],[Ap(2) Bp(2) Cp(2) Ap(2)])
% hold on
% plot(xw(:,1),xw(:,2),'*');
% 
% Area=abs(sA*(tB-tC)+sB*(tC-tA)+sC*(tA-tB))/2.0;
% Area2 = sum(xw(:,3));

s = xw(:,1); t = xw(:,2); 
w = xw(:,3);

theta = theta_x(s,t); phi = phi_x(s,t);
X = sph2cart(theta_O,phi_O,theta,phi);
X = [X(:,1) X(:,2)];
W = w.*J(s,t);

% Debug : plot final quad
% plot([A(1) B(1) C(1) A(1)],[A(2) B(2) C(2) A(2)]);
% hold on
% axis equal
% plot(X(:,1),X(:,2),'*');

end