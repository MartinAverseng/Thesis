N = 10000;
h = 1/N;
x = [0:h:1];
tic
% random coefficients
a = randi(2,N,1);
D = spdiags(a,0,N,N);
% Gradient matrix
u = ones(N,1);
G = spdiags([u/h -u/h],[0 -1],N,N-1);
% Matrix of the problem
A = G'*D*G;
% RHS
F = -G'*D*u;

% Corrector
chi = A\F;
toc
%[chi,FLAG,RELRES,ITER,RESVEC] = pcg(A,F,1e-9,5000);
[chi2] = pcg(A,F,1e-9,100);
[chi3] = pcg(A,F,1e-9,50);
plot(x,[0,chi',0],'r',x,[0,-chi2'+chi',0],'b',x,[0,-chi3'+chi',0],'g')
% Aeff
aeff = h*sum(D*G*chi + D*u)
aeff_ex = 1/(0.5*1+0.5*1/2)