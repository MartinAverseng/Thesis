

addpath('../Gypsilab/OpenMsh');
addpath('../Gypsilab/OpenMmg');
addpath('../Gypsilab/OpenDom');
addpath('../Gypsilab/OpenFem');
addpath('../Gypsilab/OpenHmx');

% Parameters
N   = 10
tol = 1e-4
typ = 'P0'
gss = 3;

% Incident wave
PW = @(X) ones(size(X,1),1);

% Mesh the disk and the half sphere S2
mesh = mshSquare(N,[1 1]);
mesh2 = mesh;
mesh2.wgt =1+30000000*rand(size(mesh.elt,1),1);

% Domain
sigma  = dom(mesh,gss);   
sigma2 = dom(mesh2, gss);

%%% PREPARE OPERATOR
disp('~~~~~~~~~~~~~ PREPARE OPERATOR ~~~~~~~~~~~~~')

% Projected Green kernel  
Gxy = @(X,Y) femGreenKernel(X,Y,'[1/r]',0);

% Finite elements on the sphere
Vh = fem(mesh,typ);

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' G(x,y) psi(y) dx dy 
tic
S  = 1/(4*pi) .* integral(sigma,sigma,Vh,Gxy,Vh);
S2 = 1/(4*pi) .* integral(sigma2,sigma2,Vh,Gxy,Vh);
toc

% Regularization
tic
Sr   = 1/(4*pi) .* projRegularize(sigma,  sigma,  Vh, '[1/r]', Vh);
Sr2  = 1/(4*pi) .* projRegularize(sigma2, sigma2, Vh, '[1/r]', Vh);
S  = S + Sr;
S2 = S2 + Sr2;
toc

N = size(S,1);
u = rand(N,1);
v = u./mesh2.wgt;

v'*S2*v 
u'*S*u

