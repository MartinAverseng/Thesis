function [theta_y,phi_y] = angles(Y,theta_x,phi_x)
% angles(Y,X) computes the angles of the vector Y = [y1,y2,y3] on the
% sphere (if not on the sphere, by default, Y will be normalized to
% Y/||Y||) in the polar basis with X as the north pole, where x has the 
% angular coordinates theta_x,phi_x in the polar basis with [0;0;1] as the
% north pole. 

normY = sqrt( Y(:,1).^2 +  Y(:,2).^2 +  Y(:,3).^2);
Y(:,1) = Y(:,1)./normY;
Y(:,2) = Y(:,2)./normY;
Y(:,3) = Y(:,3)./normY;

% Coordinates of ex, ey and ez in the basis er, etheta, ephi:
e1 = [sin(theta_x).*cos(phi_x), cos(theta_x).*cos(phi_x), -sin(phi_x)];
e2 = [sin(theta_x).*sin(phi_x), cos(theta_x).*sin(phi_x), cos(phi_x)];
e3 = [cos(theta_x), -sin(theta_x) 0];

% coordinates of y in the basis er etheta ephi:

y_er = Y(:,1).*e1(:,1) + Y(:,2).*e2(:,1) + Y(:,3).*e3(:,1);
y_etheta = Y(:,1).*e1(:,2) + Y(:,2).*e2(:,2) + Y(:,3).*e3(:,2);
y_ephi = Y(:,1).*e1(:,3) + Y(:,2).*e2(:,3) + Y(:,3).*e3(:,3);

% angles of Y in the basis (O,x):
theta_y = acos(y_er);
phi_y = atan2(y_ephi,y_etheta);







end

