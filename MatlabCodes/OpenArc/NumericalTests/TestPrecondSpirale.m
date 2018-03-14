% Test du préconditionnement sur un arc non plat.

clear all
close all
clc;

N = 1500;
x = @(s)(exp(s).*cos(5*s));
y = @(s)(exp(s).*sin(5*s));
I = [-1,1];
curve = SimpleCurve(x,y,I);
% 
% curve = polygonCurve( [0 0;0 1; 1 1],false);


meshAdapt = MeshCurve(curve,N,@cos,[-pi,0]);
L = sum(meshAdapt.length);
Vh =  weightedFEspace(meshAdapt,'P1','1/sqrt(1-t^2)',5);
M = full(Vh.Mass);

Wh =  weightedFEspace(meshAdapt,'P1','sqrt(1-t^2)',5);
MM = full(Wh.Mass);
k = 0;
X = Vh.gaussPoints;
Somega = singleLayer(k,Vh,X,{'full',true},pi);
Sgalerk = Somega.galerkine(Vh,'U');

dM = full(Wh.dMass);

dMameliore = dM + (1/log(2))^2*M*ones(size(M))/sum(sum(M))*M;
cond(4*M^(-1)*full(Sgalerk)*M^(-1)*dMameliore*M^(-1)*full(Sgalerk))


[P1,D1] = eig(dM,M);
[P2,D2] = eig(full(Sgalerk),M);
fprintf('norm( dM - M*P1*D1*P1^(-1) ) = %s \n',num2str(norm( dM - M*P1*D1*P1^(-1))));
idx = find(diag(D1) < 1);
D1(idx,idx) = 1/log(2);
sqrtdM = M*P1*sqrt(D1)*P1^(-1);
cond(2*M^(-1)*sqrtdM*M^(-1)*full(Sgalerk))
disp(cond(full(Sgalerk)));

eig(4*M^(-1)*full(Sgalerk)*M^(-1)*dMameliore*M^(-1)*full(Sgalerk))