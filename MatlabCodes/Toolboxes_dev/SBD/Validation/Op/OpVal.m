
%% Log 
close all
clear all;
clc;
% Initialize problem
M = 10^6;
N = 10^6;
[X,Y,V,Xaxis,Yaxis] = GridAndCloud(M,N);
tol = 1e-3;
% Define kernel : 
R = 10;
logK = LogKernel(R);
% Define the operator
A = Op(X,logK,Y,'tol',tol,'a_factor',1);
A = A.balanceTime;
A.show;
% Compute potential : 
q = A*V;

% Validate : 
tVal = 2; % seconds
assert(A.validate(V,q,tVal,tol)<tol);
figure
imagesc(Xaxis,Yaxis,reshape(real(q),length(Xaxis),length(Yaxis)));

%% Y0
clc;
close all;
clear all;
% Initialize problem
M = 5*10^4;
N = 10^4;
[X,Y,V,Xaxis,Yaxis] = GridAndCloud(M,N);
tol = 1e-3;
% Define kernel : 
R = 30;
hkern = H0Kernel(R);
A = Op(X,hkern,Y,'tol',tol,'a_factor',5);
q = A*V;
tVal = 2; % seconds
assert(A.validate(V,q,tVal,tol)<tol);
imagesc(Xaxis,Yaxis,reshape(abs(q),length(Xaxis),length(Yaxis)));

%% 
clear all;
close all;
clc;

Y = uniformCircle([0,0],1,50000);
X = uniformCircle([0,0],1,50000);

k = 45;
kernel  = Y0Kernel(k);
t1 = tic;
A1 = Op(Y,kernel,Y,'tol',1e-3,'rMax',2.2,'a_factor',8);
t1 = toc(t1);
disp(t1);
A1.show;
e1 = A1.validate;
t2 = tic;
A2 = A1.update_X(X);
t2 = toc(t2);
A2.show;
disp(t2);
e2 = A2.validate;


%% 
clc 
close all
disp('success');