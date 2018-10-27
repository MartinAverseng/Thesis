%% Test des résultats de préconditionnement avec la méthode de Marion Darbas. 
Main;
N = 5000;
q = 3;
r = 1; % do not change r
curve = circle(r,[0,0]);
mesh = MeshCurve(curve,N);
tol = 1e-5;

Vh = FEspace(mesh,'P1',q);
k = 25;

X =  R2toRfunc(@(Z)(Z(:,1)));
Y =  R2toRfunc(@(Z)(Z(:,2)));
THETA = atan2(Y,X);
planeWave = exp(-1i*k*X);
dr_planeWave = -1i*k*cos(THETA)*planeWave;
dr_inc = Vh.Pi_h(dr_planeWave);

S = Galerkine_singleLayer(k,Vh,Vh,'tol',tol,'a_factor',5);
g = Vh.Pi_h(planeWave);
figure
g.show;
drawnow


Mass = Vh.Mass.concretePart;
dMass = Vh.dMass.concretePart;
Np = 8;
theta = pi/3;
keps = k+1i*0.7*k^(1/3);
PrecDarbas = @(x)(-2*padePrecondDarbas(x,Np,theta,keps,Mass,dMass));

l = Vh.secondMember(planeWave);
dn_u = explicitDtN_planeWave_0inc(k,150);
phi_EXACT = Vh.Pi_h(dn_u);
figure
show(real(phi_EXACT));
hold on
lambda0 = padePrecondDarbas(l.concretePart,Np,theta,keps,Mass,dMass) - dr_inc.v;
show(real(FE_func(Vh,lambda0)))
t1 = tic;
[phi_BEM1,flag1,relres1,iter1,resvec1]  = variationalSol(S,-l,[],1e-8,300);
t1 = toc(t1);
t2 = tic;
[phi_BEM2,flag2,relres2,iter2,resvec2]  = variationalSol(S,-l,[],1e-8,300,PrecDarbas,[],lambda0);
t2 = toc(t2);
figure
semilogy(1:length(resvec1),resvec1,'-o','DisplayName',['without preconditioner : ' num2str(t1) ' s']);
hold on
semilogy(1:length(resvec2),resvec2,'-o','DisplayName',['with preconditionner : ' num2str(t2) ' s']);
legend show
drawnow; 

beep
figure
show(real(phi_BEM2));
hold on;

phi_EXACT = Vh.Pi_h(dn_u);
show(real(phi_EXACT))


