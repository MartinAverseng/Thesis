% Test 1 
a = 0;
b = 1;
rho = besselJroots(0,1);
derivative = @(x)(1./x);
C = sqrt(2)./(abs(besselj(1,rho(:))));
f = scalProds(a,b,C,derivative,rho);
assert(abs(f/C + 0.415829)<5e-6);


% Test 2
derivative = @(x)(1./x);
a =0.1;
b = 0.9;
rho = besselJroots(0,100);
C = sqrt(2)./(abs(besselj(1,rho(:))));
resolNGL = 20;
NGL = resolNGL*length(rho);
[x,w] = Gauss_Legendre1D(NGL,a,b);
ff = w.*x.*derivative(x);
JJ = diag(C)*-besselj(1,rho(:)*x(:)');
f_eiVal = JJ*ff;
f_ei = scalProds(a,b,C,derivative,rho);
assert(norm(f_ei - f_eiVal,'inf')<1e-8);

% Test 3
derivative = @(x)(bessely(1,x));
a =0.1;
b = 0.9;
rho = besselJroots(2000,0);

C = sqrt(2)./(abs(besselj(1,rho(:))));
t1 = tic;
resolNGL = 50;
NGL = resolNGL*length(rho);
[x,w] = Gauss_Legendre1D(4000,a,b);
ff = w.*x.*derivative(x);
JJ = diag(C)*-besselj(1,rho(:)*x(:)');
f_eiVal = JJ*ff;
t1 = toc(t1);
t2 = tic;
f_ei = scalProds(a,b,C,derivative,rho);
t2 = toc(t2);
assert(norm(f_ei - f_eiVal,'inf')<1e-8);
assert(t2<t1)

disp('success');


