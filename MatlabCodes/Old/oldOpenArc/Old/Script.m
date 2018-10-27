%% Create a mesh with non-uniform density and assemble the Somega

close all;
clear all;
clc;

% \int_{w(x)dx(w(x) dx u)}

N = 1000;
x = @(t)(t);
y = @(t)(0*t);
I = [-1,1];
curve = SimpleCurve(x,y,I);
meshAdapt = MeshCurve(curve,N,@cos,[-pi,0]);
L = sum(meshAdapt.length);
Vh =  weightedFEspace(meshAdapt,'P1',@(s)(1./sqrt(s.*(2 - s))),[1:2 meshAdapt.nseg-(1:-1:0)],'quadrature',5);
%Wh =  weightedFEspace(meshAdapt,'P1',@(s)(sqrt(s.*(2 - s))),[1 meshAdapt.nseg],'quadrature',5);

%K = Wh.dMass;
%H = sqrt(full(K));
k = 0;
X = Vh.gaussPoints;
%S = singleLayer(k,Wh,X,'a_factor',5);
%Sgalerk = S.galerkine(Wh,'U');
sVert = meshAdapt.sVertices;
Sw = singleLayer(k,Vh,X,'a_factor',3,'tol',1e-6);
Swgalerk = Sw.galerkine(Vh,'U');



n = 0;
Tn = R2toRfunc(@(Z)(chebyshevT(n,Z(:,1))));
%w = R2toRfunc(@(Z)(sqrt(1-Z(:,1).^2)));
lambda = Vh.secondMember(Tn);
% Test of the K matrix
%mu = K\lambda;
%figure;
%mu.show;
%hold on;
%show(Vh.Pi_h(Tn/n^2))

lambdaProj = Vh.Pi_h(Tn);
figure;
plot(X(:,1),real(Sw*lambdaProj));
hold on;
if n==0
    plot(X(:,1),ones(size(X,1),1)*log(2)/2);
else
    plot(X(:,1),1/(2*n)*Tn(X));
end

mu = variationalSol(Swgalerk,lambda,[],1e-8,100);


figure;
plot(mu);
hold on
if n==0
    plot(Vh.Pi_h(2/log(2)*Tn))
else
    plot(Vh.Pi_h(2*n*Tn))
end