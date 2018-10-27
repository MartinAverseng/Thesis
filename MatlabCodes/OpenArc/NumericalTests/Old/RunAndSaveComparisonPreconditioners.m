
a_factor = 15;
ns = [50 100 200 400 800 1600 3200];
t0_save = zeros(size(ns));
t1_save = zeros(size(ns));
t2_save = zeros(size(ns));
t3_save = zeros(size(ns));
t3bis_save = zeros(size(ns));
t4_save = zeros(size(ns));
nit0 = zeros(size(ns));
nit1 = zeros(size(ns));
nit2 = zeros(size(ns));
nit4 = zeros(size(ns));
tAssembleS_save = zeros(size(ns));
tN_save = zeros(size(ns));
i = 0;
Z0 = [-0.9,0.01];
for N = ns
    i = i+1;
    tN = tic;
    u0 = R2toRfunc(@(Z)(log(sqrt((Z(:,1)-Z0(:,1)).^2 + (Z(:,2)-Z0(:,2)).^2)).*sin(5*Z(:,1))));
    mesh = mesh.remesh(N);
    Wh = weightedFEspace(mesh,'P1','sqrt(1-t^2)',5);
    dM = Wh.dMass.concretePart;
    Vh = weightedFEspace(mesh,'P1','1/sqrt(1-t^2)',5);
    M = Vh.Mass.concretePart;
    tAssembleS = tic;
    Somega = singleLayer(0,Vh,[],{'a_factor',a_factor,'verbose',0});
    tAssembleS_save(i) = toc;
    Sgalerk = Somega.galerkine(Vh,'U');
    
    l = Somega.Vh.secondMember(u0);
    % No Prec
    t0 = tic;
    [lambda0,~,~,~,resvec0]  = variationalSol(Sgalerk,l,40,1e-8,N);
    t0_save(i) = toc(t0);
    nit0(i) = length(resvec0);
    % Prec S_0
    t1 = tic;
    Prec1 = Sgalerk.concretePart;
    [lambda1,~,~,~,resvec1]  = variationalSol(Sgalerk,l,40,1e-8,N,Prec1);
    nit1(i) = length(resvec1);
    t1_save(i) = toc(t1);
    % Prec sqrt of Delta à la volée
    t2 = tic;
    [L,U] = lu(M);
    Prec2 = @(u)(U\(L\TrefethenSqrt(4*dM,3,U\(L\u),M,4,4.5*Vh.ndof^2)));
    [lambda2,~,~,~,resvec2]  = variationalSol(Sgalerk,l,40,1e-8,N,Prec2);
    t2_save(i) = toc(t2);
    nit2(i) = length(resvec2);
    % Prec sqrt of Delta full assemble
    t3 = tic;
    Minv = M^(-1);
    Mat_Prec3 = Minv*TrefethenSqrt(4*dM,3,[],M,3,4.5*Vh.ndof^2)*Minv;
    Prec3 = @(u)(Mat_Prec3*u);
    t3bis = tic;
    [lambda3,~,~,~,resvec3]  = variationalSol(Sgalerk,l,40,1e-8,N,Prec3);
    t3bis_save(i) = toc(t3bis);
    t3_save(i) = toc(t3);
    % Prec S_\omega\Delta
    t4 = tic;
    Prec4 = @(u)(-4*(M\(Sgalerk*(M\(dM*(M\u))))));
    [lambda4,~,~,~,resvec4]  = variationalSol(Sgalerk,l,40,1e-8,N,Prec4);
    nit4(i) = length(resvec4);
    t4_save(i) = toc(t4);
    tN_save(i) = toc(tN);    
end