%% Ordre de convergence simple couche segment Laplace%% Create a mesh with non-uniform density and assemble the Somega

close all;
clear all;
clc;

x = @(t)(t);
y = @(t)(0*t);
I = [-1,1];
curve = SimpleCurve(x,y,I);
k = 0;
Nh = [10 20 40 80 160 320 640];
Xplot = linspace(-1,1,200)';
X = R2toRfunc(@(Z)(Z(:,1)));
t = -0.5;
u0 = 1/2*log(2/sqrt(1-2*t*X + t^2));
lambda = (1-t*X)/(1 - 2*t*X + t^2);
%plot(Xplot,u0([Xplot 0*Xplot]));


for i = 1:length(Nh)
    N = Nh(i);
    meshAdapt = MeshCurve(curve,N,@cos,[-pi,0]);
    L = sum(meshAdapt.length);
    Vh =  weightedFEspace(meshAdapt,'P1','1/sqrt(1-t^2)',15);
    X = Vh.gaussPoints;
    Sw = singleLayer(k,Vh,X,{'a_factor',5,'tol',1e-8});
    Swgalerk = Sw.galerkine(Vh,'U');
    u0_h = Vh.secondMember(u0);
    lambda_h = variationalSol(Swgalerk,u0_h,[],1e-10,N);
    eh = lambda_h-Vh.Pi_h(lambda);
    errH_12(i) = sqrt(real((Swgalerk*eh|eh)));
    errL2(i) = sqrt(real(eh|eh));
    close all
    plot(lambda_h);
    hold on
    plot(Vh.Pi_h(lambda))
    pause(0.1);
end

figure
loglog(1./Nh(1:i),real(errH_12(1:i)),'DisplayName','H_12')
hold on;
loglog(1./Nh(1:i),real(errL2(1:i)),'DisplayName','L2w')
legend show


