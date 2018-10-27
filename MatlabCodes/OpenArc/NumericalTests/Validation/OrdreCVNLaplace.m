%% Ordre de convergence simple couche segment Laplace
%% Create a mesh with non-uniform density and assemble the Somega

close all;
clear all;%#ok
clc;

x = @(t)(t);
y = @(t)(0*t);
I = [-1,1];
curve = SimpleCurve(x,y,I);

k = 0;

Ns = [30 60 100];
n = 2;
lambda1 = R2toRfunc.Un(n);
u01 = (n+1)/2*R2toRfunc.Un(n);

[errH_12_nonreg,errL2_nonreg] = testN_OrderCV(curve,lambda1, u01, Ns,k,'quadNum',10,'fullMatrix',true);
figure
loglogTrislope(1./Ns(:),abs(errH_12_nonreg))
hold on
loglogTrislope(1./Ns(:),errL2_nonreg)