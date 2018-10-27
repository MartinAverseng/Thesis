% Test 1 
clear all
close all
clc
a = 0;
b = 1;
rho = besselJroots(0,1);
func = @(x)(log(x));
der = @(x)(1./x);
kernel = Kernel(func,der);

f = scalProds([a,b],kernel,rho);
assert(abs(f + 2*pi*Cp(rho))<5e-6);


% Test 2
% The method with Gauss-Legendre quadrature is quicker, but number of
% quadrature points might depend on 'der' (if it has oscillations mainly). 
% So we don't take the trouble to code a better function. 

clear all
close all
clc
a =0;
b = 1;
func = @(x)(log(x));
der = @(x)(1./x);
kernel = Kernel(func,der);
rho = besselJroots(0,100);
[x,w] = Gauss_Legendre1D(fix(1.5*(b-a)*length(rho)),a,b);
ff = w.*x.*der(x);
JJ = -repmat(rho(:).*Cp(rho(:)),1,length(x)).*besselj(1,rho(:)*x(:)');
f_eiVal = 2*pi*JJ*ff(:);
f_ei = scalProds([a,b],kernel,rho);
assert(norm(f_ei - f_eiVal,'inf')<1e-8);

% Test 3
clear all
close all
clc
func = @(x)(log(x));
der = @(x)(1./x);
kernel = Kernel(func,der);
a =0.1;
b = 0.9;
rho = besselJroots(0,500);

t1 = tic;
[x,w] = Gauss_Legendre1D(fix(1.5*(b-a)*length(rho)),a,b);
ff = w.*x.*der(x);
JJ = -repmat(rho(:).*Cp(rho(:)),1,length(x)).*besselj(1,rho(:)*x(:)');
f_eiVal = 2*pi*JJ*ff;
t1 = toc(t1);
t2 = tic;
f_ei = scalProds([a,b],kernel,rho);
t2 = toc(t2);
assert(norm(f_ei - f_eiVal,'inf')<1e-8);

% Test 4
clear all
close all
clc
func = @(x)(log(x));
der = @(x)(1./x);
kernel = Kernel(func,der);
a =0.01;
b = 1;
rho = besselJroots(0,500);

t1 = tic;
[x,w] = Gauss_Legendre1D(fix(1.5*(b-a)*length(rho)),a,b);
ff = w.*x.*der(x);
JJ = -repmat(rho(:).*Cp(rho(:)),1,length(x)).*besselj(1,rho(:)*x(:)');
f_eiVal = 2*pi*JJ*ff;
t1 = toc(t1);
t2 = tic;
[f_ei,normH10,scal01] = scalProds([a,b],kernel,rho);
t2 = toc(t2);
assert(norm(f_ei - f_eiVal,'inf')<1e-8);

disp('success');





