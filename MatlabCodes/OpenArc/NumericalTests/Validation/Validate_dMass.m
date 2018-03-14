%% We test the Mass matrix. (wdx)^2
% m(u,v) = \int_{gamma} (dx(wdxu))v = - \int_{\Gamma} w dx u dx v
clear all
close all
clc;

N = 500;
x = @(t)(t);
y = @(t)(0*t);
I = [-1,1];
curve = SimpleCurve(x,y,I);
meshAdapt = MeshCurve(curve,N,@cos,-pi,0);
L = sum(meshAdapt.length);
Vh =  weightedFEspace(meshAdapt,'P1',@(s)(1./sqrt(s.*(L - s))),[1 meshAdapt.nseg],'quadrature',5);
Wh =  weightedFEspace(meshAdapt,'P1',@(s)(sqrt(s.*(L - s))),[1 meshAdapt.nseg],'quadrature',5);

dM = Wh.dMass;
n = 30;
Tn = R2toRfunc.Tn(n);

f = Vh.secondMember(-n^2*Tn);
u = variationalSol(-dM,f,[],1e-8,N);
u.show;
hold on
show(Vh.Pi_h(Tn));