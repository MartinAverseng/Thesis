function [] = plotOnSphere(theta,phi)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

sphere; 
hold on;

x = sin(theta).*cos(phi);
y = sin(theta).*sin(phi);
z = cos(theta);

plot3(x,y,z,'.','Markersize',20);

end

