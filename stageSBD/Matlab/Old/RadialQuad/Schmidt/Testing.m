% Test the BesselQuadSchmidt function
clear all;
close all;

D = 2;
c = Inf;


a = 0.008;
tt = 0;
b = 1-tt*a;
G = @(x)(log(x));
tol = 1e-6;

[ alpha,bessZs,fApprox,resL2 ] = BesselQuadSchmidt( D,c,G,a/2,b,tol );