%% Validation script for function 'radialQuad'. 


%% Test 1
close all;
clear all;
clc;
% Check if function works in the most obvious case
a = 0.1;
R = 1;
kernel = LogKernel(R);
tol = 1e-1;

out  = kernel.radialQuadKernel(a,tol);
disp(out);
show(out,1);

%% Test 2
% Check if function is stable for smaller value of a.
% It should issue a warning, since the tolerance is unreachable. 
% It also tests the feature according to which the value of 'a' is updated
% inside the function when the gram Matrix is ill-conditioned. ( a = 3/4*a)

close all;
clear all;
clc;
a = 0.01;
R = 1;
kernel = LogKernel(R);
tol = 1e-12;

out  = kernel.radialQuadKernel(a,tol);
disp(out);
show(out);

%% Test 3 (helmholtz)
close all;
clear all;
clc;
a = 0.1;
k = 15;
kernel = Y0Kernel(k);
tol = 1e-3;

out  = kernel.radialQuadKernel(a,tol);
disp(out);
show(out);

%% Test 4 (helmholtz) More coefficients
% This tests the feature of specifying a Pmax different than default value.

close all;
clear all;
clc;
a = 0.001;
k = 50;
kernel = Y0Kernel(k);
tol = 1e-3;
Pmax = 1700;

out  = kernel.radialQuadKernel(a,tol,'Pmax',Pmax);
disp(out);
show(out);

%% Test 5 (helmholtz) Big frequency, large a
% Tests non-zero start frequency feature. 
close all;
clear all;
clc;
a = 0.5;
k = 500;
R = nextY0root(k);
kernel = Y0Kernel(R);
tol = 1e-3;

out  = kernel.radialQuadKernel(a,tol);
disp(out);
show(out);


%% Test 6 (Non-homogenous kernel without multi-Dirichlet conditions)

close all;
clear all;
clc;
a = 0.01;
func = @(x)(sin(10*x));
derivative = @(x)(10*cos(10*x));
kernel = Kernel(func,derivative);
tol = 1e-5;

out  = kernel.radialQuadKernel(a,tol);
disp(out);
show(out);

%% Test 7 : Monitor derivative of target kernel

close all;
clear all;
clc;


a = 0.01/1.05;
kernel = LogKernel(1);
tol = 1e-5;

out  = kernel.radialQuadKernel(a,tol,'monitorDerivative',true);
disp(out);
showDer(out);

kernel = Kernel(@(x)(1./x.^2),@(x)(-2./x.^3));
out  = kernel.radialQuadKernel(a,tol);
disp(out);
show(out);





%%
clc
close all
disp('success')