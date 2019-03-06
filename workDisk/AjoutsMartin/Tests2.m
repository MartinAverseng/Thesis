% Changement de coordonnées sphériques:
close all
clear all;
clc;

thetax = pi/2*rand(1,1);
phix = 2*pi*rand(1,1);
% thetax et phix dans la base polaire de [0;0;1]

thetay = pi*rand(10,1);
phiy = 2*pi*rand(10,1);
% thetay et phiy dans la base polaire de x. 

erethetaephix = [[sin(thetax).*cos(phix);sin(thetax).* sin(phix);cos(thetax)]...
    [cos(thetax).*cos(phix); cos(thetax).*sin(phix); -sin(thetax)]...
    [-sin(phix); cos(phix); 0]];


B = erethetaephix';

erx = erethetaephix(:,1);
ethetax = erethetaephix(:,2);
ephix = erethetaephix(:,3);
x = erx;

ery = [sin(thetay).*cos(phiy) sin(thetay).* sin(phiy) cos(thetay)];
ethetay = [cos(thetay).*cos(phiy) cos(thetay).*sin(phiy) -sin(thetay)];
ephiy = [-sin(phiy); cos(phiy); 0*phiy];

% On obtient les coordonnées cartésiennes de y:
[y] = sph2cart(thetax,phix,thetay,phiy);
% On exprime maintenant les angles de y dans la base (O,z):
y1 = y(:,1);
y2 = y(:,2);
y3 = y(:,3);
thetay2 = acos(y3);
phiy2 = atan2(y2,y1);

% Autre manière de calculer les coords de y:

ery2 = [sin(thetay2).*cos(phiy2) sin(thetay2).* sin(phiy2) cos(thetay2)];



sphere
hold on
vectarrow(x,x + erx)
hold on
vectarrow(y(7,:),y(7,:) + y(7,:))
hold on
vectarrow(x,x + ethetax)
hold on
vectarrow(x,x + ephix)
hold on
vectarrow( y(7,:)+  y(7,:),y(7,:)+  y(7,:)+ery2(7,:));


axis equal



