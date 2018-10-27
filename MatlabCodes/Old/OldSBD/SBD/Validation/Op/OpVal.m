
%% Log 
close all
clear all;
clc;
% Initialize problem
M = 2*10^5;
N = 10^2;
[X,Y,V,Xaxis,Yaxis] = GridAndCloud(M,N);
tol = 1e-3;
% Define kernel : 
R = 10;
y0kern = LogKernel(R);
% Choose 'a' paramter
Ngeom = sqrt(N*M);
a = 1/(1.5*sqrt(Ngeom)); 
% Optimized empirically for random clouds
b = 1; 
% Define the operator
A = Op(X,Y,y0kern,a,b,tol);

% Compute potential : 
q = A*V;

% Validate : 
tVal = 5; % seconds
assert(A.validate(V,q,tVal,tol)<tol);
figure
imagesc(Xaxis,Yaxis,reshape(real(q),length(Xaxis),length(Yaxis)));

%% Y0
clc;
close all;
clear all;
% Initialize problem
M = 5*10^5;
N = 10^5;
[X,Y,V,Xaxis,Yaxis] = GridAndCloud(M,N);
tol = 1e-3;
% Define kernel : 
R = 30;
y0kern = Y0Kernel(R);
j0kern = J0Kernel(R);
% Choose 'a' paramter
Ngeom = sqrt(N*M);
a = 1/(sqrt(Ngeom)); 
% Optimized empirically for random clouds
b = 1; 
% Define the operator
A1 = Op(X,Y,y0kern,a,b,tol);
A2 = Op(X,Y,j0kern,0,b,tol);
% Compute potential : 
q1 = A1*V;
tVal = 5; % seconds
assert(A1.validate(V,q1,tVal,tol)<tol);
q2 = A2*V;
tVal = 5; % seconds
assert(A2.validate(V,q2,tVal,tol)<tol);

q = q1 - 1i*q2;
figure
imagesc(Xaxis,Yaxis,reshape(abs(q),length(Xaxis),length(Yaxis)));


%% 
clc 
close all
disp('success');