run('Main.m')
% Adds all functions of this folder and subfolders into  the path. 

%% Example 1 : create an operator for the Laplace kernel : 
% We want to compute q = A*V where V is the potential and A(i,j) =
% log(R|x_i - y_j|) where x and y are two random cloud of points, R is some
% constant. 
clear all;
close all;
clc;

% Generating random data
Nx = 10^5;
Ny = 10^5;
[X,Y,V] = randomCloud(Nx,Ny);
tol = 1e-3;

% Define kernel
R = 10; 
kernel = LogKernel(R);

% Choose a value for the parameter a. 
Ngeom = sqrt(Nx*Ny);
a = 2/(sqrt(Ngeom)); % This value may be chosen by user, or optimized 
% numerically by a dichotomy strategy (balancing close and far field
% constants)
% Here the choice in 1/sqrt(N) is dictated by distribution density 

b = 1; % For usual kernels (Log, J0 and Y0), always set b to 1. 
% For user-defined kernel, it can be worth trying to choose b <1 although
% no theory available to justify convergence...


% Build operator
A = Op(X,Y,kernel,a,b,tol);
B = Op(X,X,kernel,a,b,tol);
A.show;

% Perform matrix-vector product and validate
q = A*V;
tVal = 5; % Number of seconds allowed for computing exact full product on a sample of the points X. 
err = A.validate(V,q,tVal,tol); % this line issues a warning if the tolerance is not fullfilled.  
% variable err contains the maximal error found during tVal between the
% approximated value and the exact solution. 


%% Example 2 : create an operator for the Helmholtz kernel : 
% We want to compute q = A*V where V is the potential and A(i,j) =
% Y_0(R|x_i - y_j|) - i*J0(R|x(i) - y(j)|) where x and y are two random cloud of points, R is some
% constant. 
clear all;
close all;
clc;

% Generating random data
Nx = 10^4;
Ny = 10^4;
[X,Y,V] = randomCloud(Nx,Ny);
tol = 1e-3;

% Define kernel
R = 10; 
y0kern = Y0Kernel(R);
j0kern = J0Kernel(R);
% Choose a value for the parameter a. 
Ngeom = sqrt(Nx*Ny);
a = 2/(sqrt(Ngeom)); % This value may be chosen by user, or optimized 
% numerically by a dichotomy strategy (balancing close and far field
% constants)
% Here the choice in 1/sqrt(N) is dictated by distribution density 

b = 1; % For usual kernels (Log, J0 and Y0), always set b to 1. 
% For user-defined kernel, it can be worth trying to choose b <1 although
% no theory available to justify convergence...


% Build operator
A1 = Op(X,Y,y0kern,a,b,tol);
A2 = Op(X,Y,j0kern,0,1,tol);
A1.show;

% Perform matrix-vector product and validate
q1 = A1*V;
q2 = A2*V;
q = q1 - 1i*q2;
tVal = 5; % Number of seconds allowed for computing exact full product on a sample of the points X. 
err1 = A1.validate(V,q1,tVal,tol); % this line issues a warning if the tolerance is not fullfilled.  
err2 = A2.validate(V,q2,tVal,tol); 
% variable err contains the maximal error found during tVal between the
% approximated value and the exact solution. 

%% Example 3 : create an operator for an arbitrary kernel 
% We want to compute q = A*V where V is the potential and A(i,j) =
% G(|x_i - y_j|) where x and y are two random cloud of points, G is some
% kernel. 
clear all;
close all;
clc;

% Generating random data
Nx = 10^5;
Ny = 10^5;
[X,Y,V] = randomCloud(Nx,Ny);
tol = 1e-3;

% Define kernel
%kernel = HelmholtzPerturb(100);
kernel = Kernel(@(X)(cos(X) - (-cos(1) - sin(1))/4*X.^2),@(X)(-sin(X) - (-cos(1) - sin(1))/2*X));
% Choose a value for the parameter a. 
Ngeom = sqrt(Nx*Ny);
a = 2/(sqrt(Ngeom)); % This value may be chosen by user, or optimized 
% numerically by a dichotomy strategy (balancing close and far field
% constants)
% Here the choice in 1/sqrt(N) is dictated by distribution density 

b = 1; % For usual kernels (Log, J0 and Y0), always set b to 1. 
% For user-defined kernel, it can be worth trying to choose b <1 although
% no theory available to justify convergence...


% Build operator
A = Op(X,Y,kernel,a,b,tol);
A.show;

% Perform matrix-vector product and validate
q = A*V;
tVal = 5; % Number of seconds allowed for computing exact full product on a sample of the points X. 
err = A.validate(V,q,tVal,tol); % this line issues a warning if the tolerance is not fullfilled.  
% variable err contains the maximal error found during tVal between the
% approximated value and the exact solution. 


%% Example 4 
% Computes a beautiful figure looking like a mandala
clear all;
close all;
clc;

% Generating data
Nx = 5*10^6;
Ny = 10;
Ny2 = 50;
[X,Y,V,Xaxis,Yaxis] = GridAndUniformCircleData(Nx,Ny);
[X2,Y2,V2,Xaxis2,Yaxis2] = GridAndUniformCircleData(Nx,Ny2);
Y = 1.5*Y;
Y2 = 0.5*Y2;
tol = 1e-3;

% Define kernel
R = 50; 
y0kern = Y0Kernel(R);
j0kern = J0Kernel(R);

R2 = 200; 
y0kern2 = Y0Kernel(R2);
j0kern2 = J0Kernel(R2);

% Choose a value for the parameter a. 
Ngeom = sqrt(Nx*Ny);
a = 1/(1.5*sqrt(Ngeom))*5; % This value may be chosen by user, or optimized 
% numerically by a dichotomy strategy (balancing close and far field
% constants)
% Here the choice in 1/sqrt(N) is dictated by distribution density 

b = 1; % For usual kernels (Log, J0 and Y0), always set b to 1. 
% For user-defined kernel, it can be worth trying to choose b <1 although
% no theory available to justify convergence...


% Build operator
A1 = Op(X,Y,y0kern,a,b,tol);
A2 = Op(X,Y,j0kern,a,b,tol);
A3 = Op(X2,Y2,y0kern2,a,b,tol);
A4 = Op(X2,Y2,j0kern2,a,b,tol);
A1.show;

% Perform matrix-vector product and validate
q1 = A1*V;
q2 = A2*V;
q3 = A3*V2;
q4 = A4*V2;
q = q1 - 1i*q2;
q(abs(q) > 0.8)=0.8;
r = q3 - 1i*q4;
r(abs(r) > 0.8)=0.8;
tVal = 5; % Number of seconds allowed for computing exact full product on a sample of the points X. 
err1 = A1.validate(V,q1,tVal,tol); % this line issues a warning if the tolerance is not fullfilled.  
err2 = A2.validate(V,q2,tVal,tol); 
% variable err contains the maximal error found during tVal between the
% approximated value and the exact solution. 

figure
imagesc(Xaxis,Yaxis,reshape(abs(q + r),length(Xaxis),length(Yaxis)));

%% Example 5
% Far field Helmholtz

clear all;
close all;
clc;

% Generating random data
Nx = 10^6;
Ny = 10^6;
[X,Y,V] = TwoFarClouds(Nx,Ny);
tol = 1e-3;

% Define kernel
R = 50; 
y0kern = Y0Kernel(R);
j0kern = J0Kernel(R);
% Choose a value for the parameter a. 
Ngeom = sqrt(Nx*Ny);
a = 0.5; % Big value

b = 1; % For usual kernels (Log, J0 and Y0), always set b to 1. 
% For user-defined kernel, it can be worth trying to choose b <1 although
% no theory available to justify convergence...


% Build operator
A1 = Op(X,Y,y0kern,a,b,tol,'noCloseField',true); % we know that no close interactions are required
A2 = Op(X,Y,j0kern,0,b,tol,'noCloseField',true); % we know that no close interactions are required
A1.show;

% Perform matrix-vector product and validate
q1 = A1*V;
q2 = A2*V;
q = q1 - 1i*q2;
tVal = 5; % Number of seconds allowed for computing exact full product on a sample of the points X. 
err1 = A1.validate(V,q1,tVal,tol); % this line issues a warning if the tolerance is not fullfilled.  
err2 = A2.validate(V,q2,tVal,tol); 
% variable err contains the maximal error found during tVal between the
% approximated value and the exact solution. 







