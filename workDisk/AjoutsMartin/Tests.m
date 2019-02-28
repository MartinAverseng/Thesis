% Test the integrations:

clear all
close all
clc;

% Reference triangle

addpath('../Gypsilab/OpenMsh');
addpath('../Gypsilab/OpenMmg');
addpath('../Gypsilab/OpenDom');
addpath('../Gypsilab/OpenFem');
addpath('../Gypsilab/OpenHmx');

Aref = [0 0]; Bref = [1 0]; Cref = [0 1];
% Integration of a polynomial
P = @(x,y)(x.^3.*y.^4 + x.^2.*y.^1);
Ngauss = 3;
err = abs(gaussQuadTri(P,Aref,Bref,Cref,Ngauss) - integral2tri(P,Aref,Bref,Cref));
disp(err)
% Integration of another function
f = @(x,y)(sin(x).*cos(y));
err = abs(gaussQuadTri(f,Aref,Bref,Cref,Ngauss) - integral2tri(f,Aref,Bref,Cref));
disp(err)

% Some other triangle
A = randn(2,1); B = randn(2,1); C = randn(2,1);
% Integration of a polynomial
P = @(x,y)(x.^3.*y.^4 + x.^2.*y.^1);
err = abs(gaussQuadTri(P,A,B,C,Ngauss) - integral2tri(P,A,B,C));
disp(err)
% Integration of another function
err = abs(gaussQuadTri(f,A,B,C,Ngauss) - integral2tri(f,A,B,C));
disp(err)

% Explicit integral of a smooth function and validation
Aref = [0 0]; Bref = [1 0]; Cref = [0 1];
Xfar = [10 10]; % Far away
freg = @(x,y)(1./sqrt((Xfar(1) - x).^2 + (Xfar(2) - y).^2));
err = abs(exactIntRm1Tri(Aref,Bref,Cref,Xfar) - integral2tri(freg,Aref,Bref,Cref));
disp(err)


disp('1/|X - Y| (we have the explicit value)')
% Singularity : on a vertex

Aref = [0 0]; Bref = [1 0]; Cref = [0 1];
Xclose = Aref;
freg = @(x,y)(1./sqrt((Xclose(1) - x).^2 + (Xclose(2) - y).^2));
err1 = abs(exactIntRm1Tri(Aref,Bref,Cref,Xclose) - integral2tri(freg,Aref,Bref,Cref));
err2 = abs(exactIntRm1Tri(Aref,Bref,Cref,Xclose) - gaussQuadTri(freg,Aref,Bref,Cref,Ngauss));
fprintf('Vertex : matlab err = %s, quad err = %s\n\n',num2str(err1),num2str(err2))


% Singularity : on an edge

Aref = [0 0]; Bref = [1 0]; Cref = [0 1];
Xclose = 0.5*(Aref + Bref);
freg = @(x,y)(1./sqrt((Xclose(1) - x).^2 + (Xclose(2) - y).^2));
err1 = abs(exactIntRm1Tri(Aref,Bref,Cref,Xclose) - integral2tri(freg,Aref,Bref,Cref));
err2 = abs(exactIntRm1Tri(Aref,Bref,Cref,Xclose) - gaussQuadTri(freg,Aref,Bref,Cref,Ngauss));
fprintf('Edge : matlab err = %s, quad err = %s\n\n',num2str(err1),num2str(err2))

% Singularity : inside the triangle

Aref = [0 0]; Bref = [1 0]; Cref = [0 1];
Xclose = 0.33*(Aref + Bref + Cref);
freg = @(x,y)(1./sqrt((Xclose(1) - x).^2 + (Xclose(2) - y).^2));
err1 = abs(exactIntRm1Tri(Aref,Bref,Cref,Xclose) - integral2tri(freg,Aref,Bref,Cref));
err2 = abs(exactIntRm1Tri(Aref,Bref,Cref,Xclose) - gaussQuadTri(freg,Aref,Bref,Cref,Ngauss));
fprintf('Inside : matlab err = %s, quad err = %s\n\n',num2str(err1),num2str(err2))

%% Special function


disp('1/|F(X) - F(Y)| (this time we trust matlab)')

Aref = [0 0]; Bref = [1 0]; Cref = [0 1];
Xclose = .5*(Aref + Bref);
x1 = Xclose(1); x2 = Xclose(2);
F1 = @(x,y)(y + cos(x));
F2 = @(x,y)(x + cos(y));
F = @(x,y)([F1(x,y) F2(x,y)]);
dFX = [-sin(x1) 1; 1 -sin(x2)];

fsing = @(x,y)(1./sqrt((F1(x1,x2) - F1(x,y)).^2 + (F2(x1,x2) - F2(x,y)).^2));


err1 = abs(gaussQuadTri(fsing,Aref,Bref,Cref,Ngauss)-integral2tri(fsing,Aref,Bref,Cref));
err2 = abs(singIntSpecial(Aref,Bref,Cref,Xclose,F1,F2,dFX,Ngauss)-integral2tri(fsing,Aref,Bref,Cref));

fprintf('edge : (brutal) quad err = %s, special method err = %s \n\n',num2str(err1),num2str(err2));


