%% Validation script for function 'radialQuad'. 


%% Test 1
close all;
clear all;
clc;
% Check if function works in the most obvious case
a = 0.1;
b = 1;
R = 1;
kernel = LogKernel(R);
tol = 1e-1;

out  = kernel.radialQuadKernel(a,b,tol);
disp(out);
show(out);

%% Test 2
% Check if function is stable for smaller value of a and using quadratures
% instead of explicit forms
close all;
clear all;
clc;
a = 0.01;
b = 1;
R = 1;
kernel = LogKernel(R);
tol = 1e-3;

out  = kernel.radialQuadKernel(a,b,tol);
disp(out);
show(out);

%% Test 2 bis
% Test if Pmax < Inf is properly handled when conflict with tolerance
% A warning should be displayed about tolerance and Pmax. 
% The error shown in the figure should be greater than the tolerance
% (shown in dashed lines)
close all;
clear all;
clc;
a = 0.001;
b = 1;
R = 1;
kernel = LogKernel(R);
tol = 1e-6;
Pmax = 3000;

out  = kernel.radialQuadKernel(a,b,tol);
disp(out);
show(out);

%% Test 3 (helmholtz)
close all;
clear all;
clc;
a = 0.1;
b = 1;
k = 15;
R = nextY0root(k);
kernel = Y0Kernel(R);
tol = 1e-3;

out  = kernel.radialQuadKernel(a,b,tol);
disp(out);
show(out);

%% Test 4 (helmholtz) More coefficients
close all;
clear all;
clc;
a = 0.001;
b = 1;
k = 50;
R = nextY0root(k);
kernel = Y0Kernel(R);
tol = 1e-6;

out  = kernel.radialQuadKernel(a,b,tol);
disp(out);
show(out);

%% Test 5 (helmholtz) Big frequency, large a

close all;
clear all;
clc;
a = 0.8;
b = 1;
k = 2000;
R = nextY0root(k);
kernel = Y0Kernel(R);
tol = 1e-3;

out  = kernel.radialQuadKernel(a,b,tol);
disp(out);
show(out);

%% Test 5 bis (helmholtz) Big frequency, reasonnable a
% This example should not display any warning about a, and the number of coefficients should be very very low 
close all;
clear all;
clc;
a = 0.5;
b = 1;
k = 2000;
R = nextY0root(k);
kernel = Y0Kernel(R);
tol = 1e-3;

out  = kernel.radialQuadKernel(a,b,tol);
disp(out);
show(out);


%% Test 6 (Very explosive function)
% Issues a lot of warning because the function is not H10. 
close all;
clear all;
clc;
a = 0.01;
b = 1;
func = @(x)(1./(x.^2) + 1./(1+x));
derivative = @(x)(-2./(x.^3) - 1./(1+x).^2);
kernel = Kernel(func,derivative);
tol = 1e-3;

out  = kernel.radialQuadKernel(a,b,tol);
disp(out);
show(out);

%% 
clc
close all
disp('success')