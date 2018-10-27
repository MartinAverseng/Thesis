%% Solving the Dirichlet problem for the Laplacian outside a open arc. 
% Our aim is to solve the boundary value problem $\Delta u = 0$ with given
% Dirichlet data on an open smooth curve in $R^2$. The solution has a singularity
% near the tips of the arc causing difficulties for numerical resolution.
% We show a method of resolution using a weighted version of the single
% layer potential and a preconditioner for the system based on an analysis
% when the arc is the open segment $\Gamma = (-1,1)$. This preconditioner
% involves the square root of a weighted version of the Laplace operator. 

%% 1°) The singularity
% We first create the curve, which is the $[-1,1]$ segment.
% _segment_ is an object of type "simpleCurve", containing the parametric
% representation of the curve

segment = unitSegment;
figure
plot(segment)
title('The unit segment')
xlabel('t');

%%
% We define a mesh on the curve.
%%

Nmesh = 10;
mesh = MeshCurve(segment,Nmesh);
figure
plot(mesh);
title(sprintf('The mesh with %s segments',num2str(Nmesh)));
xlabel('t');

%%
% We can define 'P0' Finite elements on the mesh
%%

Vh = FEspace(mesh,'P0');
figure
plot(Vh.cell);
title('The ''P0'' finite element functions on the reference element');
xlabel('x');

%%
% or also 'P1' Finite elements on the mesh
%%

Vh = FEspace(mesh,'P1');
figure
plot(Vh.cell);
title('The ''P1'' finite element functions on the reference element');
xlabel('x');

%%
% Say we want to solve the Dirichlet scattering problem
% That is, let an incident field $u_i$ with $\Delta u_i = 0$ in all $\bf{R}^2$. Let us
% choose an incident wave that is constant.

u_i = R2toRfunc(1);

%%
% The scattered field $u$ satisfies $\Delta u = 0$ with $u = -u_i$ prescribed
% on both sides of the segment. Then, the scattered field is equal to
% $S\lambda$ where $\lambda$ is the jump of normal derivatives of $u$ across
% $\Gamma$ and $S$ is the single layer potential.
% On the segment, this translates to $S\lambda = - u_i$. Let us assemble the
% single layer potential. Here, the variable _S_ is an object that contains
% a method for the matrix product evaluation on each vector, without
% storing the actual coefficients of the matrix. The linear system is
% solved using Krylov method (GMRES).
% For a function $u$ defined on $\Gamma$, the secondMember routine returns
% the linear form l with coefficients
%
% $$l_i = \int_{\Gamma} u(x) \phi_i(x).$$
%
% Those coefficients are obtained via Gaussian quadratures.

k = 0; % Frequency.
S = singleLayer(k,Vh,[],{'a_factor',5}); % By default, we set the evalutation points of the
% Single layer on the gauss points of Vh.
Sgalerkine = S.galerkine(Vh,'U'); % A bilinear form with a matrix representation
l = Vh.secondMember(-u_i);

%%

lambda = variationalSol(Sgalerkine,l); % Here the gmres takes place.
figure
plot(lambda);
title('The approximated solution of S lambda = - u_i')

%%
% The exact solution of this problem is known.
t = R2toRfunc(@(Z)(Z(:,1)));
lambda_theo = -2/log(2)*1/sqrt(1-t^2);
lambda_theo_h = Vh.Pi_h(lambda_theo);
hold on;
plot(lambda_theo_h);
legend({'approximated solution','Theoretical solution'});

%%
% If we refine the mesh, the approximation converges to the actual
% solution, however quite slowly. 


S = S.remesh(300);
Sgalerkine = S.galerkine(S.Vh,'U');
l = S.Vh.secondMember(-u_i);
lambda = variationalSol(Sgalerkine,l,[],1e-8,300); % Here the gmres takes place.
figure
plot(lambda);
title('The approximated solution of S lambda = - u_i')
lambda_theo_h = S.Vh.Pi_h(lambda_theo);
hold on;
plot(lambda_theo_h);
legend({'approximated solution','Theoretical solution'});


%% Variable change, weighted operator
% The solution has a singularity that can be removed by a weighted version
% of the single layer.
% First we remesh, using a cosinusoidal repartition the mesh nodes. This is
% needed to ensure the stability of the new galerkine scheme.

repartition = @cos;
bounds = [-pi,0];
mesh = MeshCurve(segment,Nmesh,repartition,bounds);
figure;
plot(mesh);
title(sprintf('Adapted mesh with %s segments',num2str(Nmesh)));
xlabel('t');

%%
% We can now create the weighted FE space with a weight, using the
% following :

weightId = '1/sqrt(1-t^2)';
quadNum = 5;
Vh = weightedFEspace(mesh,'P1',weightId,quadNum);
figure;
Vh.showWeight;

%%
% For the weighted FE space with weight $\omega$ the mass matrix is given
% by :
%
% $$ Mij = \int_{\Gamma} \omega(x) \phi_i(x) \phi_j(x).$$
%
% Those integrals are again obtained by Gaussian quadratures, but this
% time, the nodes and weights are suited to integration against $\omega$,
% that is, a quadrature rule is computed such that
%
% $$\int_{\Gamma} \omega(x) f(x) \approx \sum_{q=1}^Q \omega_i f(x_i)$$
%
% for each segment in the mesh.


%%
% We assemble the weighted single layer :
Somega = singleLayer(k,Vh); % Same syntax.
Sgalerkine = Somega.galerkine(Vh,'U');

%%
% We now expect the solution to be a constant. To obtain the solution, we
% use the same syntax as previously.

l = Vh.secondMember(-u_i);
lambda = variationalSol(Sgalerkine,l);
figure;
plot(real(lambda),'DisplayName','Approx solution');
hold on;
lambda_theo = R2toRfunc(-2/log(2));
lambda_theo_h = Vh.Pi_h(lambda_theo);
plot(lambda_theo_h,'DisplayName','exact solution');
title('Solution of the weighted problem');
legend show;

%%
% We can solve the same problem with a finer mesh. We must also reassemble
% the single layer :
% One can observe there are still some imprecisions at the edge but this is
% actually only happening for the constant functions, and comes from a
% difficulty in the evaluation of singular integrals. This problem should
% be solvable.

Nmesh = 300;
Vh = Vh.remesh(Nmesh);
AopOpt = {'a_factor',7};
% Those options are used for the creation of the SBD operator Aop a_factor
% is a parameter to control the amount of interactions that are considered
% "close". Here, we use it to reduce the time of a matrix vector product.
Somega = singleLayer(k,Vh,[],AopOpt);
Sgalerkine = Somega.galerkine(Vh,'U');
l = Vh.secondMember(-u_i);
lambda = variationalSol(Sgalerkine,l,[],1e-8,Nmesh);
figure;
plot(real(lambda),'DisplayName','approx solution');
hold on
lambda_theo = R2toRfunc(-2/log(2));
lambda_theo_h = Vh.Pi_h(lambda_theo);
plot(lambda_theo_h,'DisplayName','exact solution');
title('Solution of the weighted problem');
legend show;

%%
% We can also solve the Dirichlet problem for more complex second-members
% We define a second member by a (finite) series of chebichev polynomials
% for which we know the exact solution. Let
%
% $$u_0(x) = c_0 + \sum_{n= 1}^N u_n T_n(x).$$
%
% Then, we have
%
% $$\lambda(x) = \frac{2}{\log(2)}c_0 + \sum_{n=1}^N 2n u_n T_n(x).$$
%

u0_series = @(n)(1./(n+1).^2 + (-1).^n./(n+1).^2);
c_u0 = 1;
lambdaSeries = @(n)(2*n.*u0_series(n));
c_lambda = 2/log(2)*c_u0;
Nmax_series = 20;
u0 = R2toRfunc(@(Z)(fromChebSeries(u0_series,Nmax_series,Z(:,1)))) + c_u0;
lambda_theo = R2toRfunc(@(Z)(fromChebSeries(lambdaSeries,Nmax_series,Z(:,1)))) + c_lambda;

lambda = variationalSol(Sgalerkine,Vh.secondMember(u0),[],1e-8,Nmesh);
figure;
plot(Vh.Pi_h(u0),'DisplayName','Second member');
hold on
plot(Vh.Pi_h(lambda_theo),'DisplayName','Exact solution solution');
plot(lambda,'LineStyle','--','DisplayName','Approximated solution');
legend show;

%% Preconditioning the linear system
%
% <html><h3>Using S0 as a preconditionner</h3></html>
%
% As we see here, the resolution of this problem requires a lot gmres
% iterations, due to the bad conditionning of the linear system to be
% solved. We could first try to precondition this system using the
% SBD decomposition of the single layer operator : S = S0 + R where R is a
% smoothing operator and S0 has a sparse matrix representation, which can
% be inverted easily. To approximate the Galerkine matrix of the opertor
% $S_0^{-1}S$, we use that if $S_0^{-1}S \lambda = u$, then, $S\lambda = S_0
% u$ so applying a test function on each side, we have :
%
% $$\int_{\Gamma} \omega(x)(S\lambda)(x)v(x) = \int_{\Gamma}\omega(x)(S_0u)(x)v(x).$$
%
% so that a consistent approximation of $u = S_0^{-1}S \lambda$ is just
% $[S_0]^{-1}[S] \Lambda$ where $\Lambda$ is the vector of coordinates of
% the projection of $\lambda$ on the mesh.

Prec = Sgalerkine.concretePart;
figure
spy(Prec);
title('Sparsity structure of the principal part of Sgalerkine')
t1 = tic;
[~,~,relres1,iter1,resvec1] = variationalSol(Sgalerkine,l,[],1e-8,Nmesh);
t1 = toc(t1);
t2 = tic;
[~,~,relres2,iter2,resvec2] = variationalSol(Sgalerkine,l,[],1e-8,Nmesh,Prec);
t2 = toc(t2);
fprintf('Solve without preconditionner : %s interations, and %s seconds \n',num2str(iter1(2)),num2str(t1));
fprintf('Solve with preconditionner : %s interations, and %s seconds \n',num2str(iter2(2)),num2str(t2));


figure
semilogy(1:length(resvec1),resvec1,'DisplayName','Without preconditioner');
hold on
semilogy(1:length(resvec2),resvec2*norm(full(Prec)),'DisplayName','With preconditionner');
legend show;
title('GMRES relative residual history');

%%
% Preconditionning by $S_0$ improves the number of iterations and the time.
% However, the number of iterations increases in both cas as the mesh is
% refined.

niter = 0;

ns = [30 50 80 100 150 300 500 1000 1500];
ngmres1 = zeros(length(ns),1);
ngmres2 = zeros(length(ns),1);
t1 = zeros(length(ns),1);
t2 = zeros(length(ns),1);
for n = ns
    niter = niter + 1;
    Somega = Somega.remesh(n);
    Sgalerkine = Somega.galerkine(Somega.Vh,'U');
    l = Somega.Vh.secondMember(-u_i);
    Prec = Sgalerkine.concretePart;
    t11 = tic;
    [~,~,~,~,resvec1] = variationalSol(Sgalerkine,l,[],1e-8,n);
    ngmres1(niter) = length(resvec1);
    t1(niter) = toc(t11);
    t22 = tic;
    [~,~,~,~,resvec2] = variationalSol(Sgalerkine,l,[],1e-8,n,Prec);
    t2(niter) = toc(t22);
    ngmres2(niter) = length(resvec2);
end

figure
subplot(1,2,1)
plot(ns,ngmres1,'DisplayName','Without preconditionner');
hold on
plot(ns,ngmres2,'DisplayName','With preconditionner');
xlabel('problem size')
ylabel('niter');
title('Iteration count')
legend show
subplot(1,2,2)
plot(ns,t1,'DisplayName','Without preconditionner');
hold on
plot(ns,t2,'DisplayName','With preconditionner');
xlabel('problem size')
ylabel('time(s)');
title('Time')
legend show




