%% Evolution du temps de calcul en fonction de N

close all;
clear all;


%% Problem parameters :

tol = 1e-4;
r = 1;
ftpy = 0;
Nvec = [1000, 5000, 10000, 20000, 30000, 40000, 50000];
for Ncharges = Nvec;
    ftpy = ftpy + 1;
    
    %% Charges
    chargesX = zeros(Ncharges,1);
    chargesY = zeros(Ncharges,1);
    NchargesTemp = 0;
    while NchargesTemp <Ncharges
        newChargeX = (rand(1,1)-0.5)*2*r;
        newChargeY = (rand(1,1)-0.5)*2*r;
        if newChargeX^2 + newChargeY^2 <= r^2
            NchargesTemp = NchargesTemp+1;
            chargesX(NchargesTemp) = newChargeX;
            chargesY(NchargesTemp) = newChargeY;
        end
    end
    chargesQ = randn(Ncharges,1)*1/Ncharges;
    
    %% Choose the parameter a
Rmax = 2;
a = sqrt(sqrt(log(1/tol)))/sqrt(Ncharges)/2.3; % To be modified later
Rmin = a*Rmax;



%% Compute the radial quadrature

D = 2;
c = Inf;
G = @(x)(log(x));
tolQuadRad = tol*a;

[ alpha,bessZs,fApprox,resL2,~,~,~,timeQuad(ftpy) ] = BesselQuadSchmidt( D,c,G,a,1,tolQuadRad,0);
close all;
P = length(alpha);


%% Circular quadrature

asympt = log(8*norm(alpha,1)/tol);
ksi = [];
w = [];
for p = 1:P
    Np = fix((exp(1)*2*bessZs(p) + asympt))+1;
    err(p) = 3*(exp(1)*bessZs(p)/Np)^Np;
    ksi = [ksi bessZs(p)*exp(1i*(0:Np-1)/Np*2*pi)];
    w = [w alpha(p)/Np*ones(1,Np)];
end
    


%% Produit matrice vecteur lointain :
tic;


% 1. Passage en fourier non uniforme

nj = Ncharges;
xj = chargesX/Rmax;
yj = chargesY/Rmax;
cj = chargesQ;
iflag = -1;

nk = length(ksi);
sk = real(ksi)';
tk = imag(ksi)';

aux = nufft2d3(nj,xj,yj,cj,iflag,1e-15,nk,sk,tk);

% 2. Retour en espace
nj = nk;
xj = sk;
yj = tk;
cj = aux.*w';
iflag = 1;

nk = Ncharges;
sk = chargesX/Rmax;
tk = chargesY/Rmax;

field = nufft2d3(nj,xj,yj,cj,iflag,1e-15,nk,sk,tk)+ G(Rmax)*sum(chargesQ);

tNuFFT(ftpy) = toc;



%% Produit proche
tic
[I,rxy] = rangesearch([chargesX(:) chargesY(:)],[chargesX(:) chargesY(:)],1.2*Rmin);
jdx = cell2mat(I')';
rxy = cell2mat(rxy')';
idx = zeros(size(jdx));
j = 1;
for i=1:length(I);
    idx(j:j+length(I{i})-1) = i;
    j = j + length(I{i});
end

% Exact value of the green kernel
Gxy = G(rxy);
Gxy(rxy < 1e-15) = 0;
% SCSD values
SCSDxy = fApprox(rxy/Rmax) + G(Rmax);

% Matrix value
val = Gxy - SCSDxy;

% Compute correction matrix
M = sparse(idx,jdx,val,Ncharges,Ncharges);

field = field + M*chargesQ;
tClose(ftpy) = toc;

end
save('timeQuad')
save('tNuFFT');
save('tClose');
load('timeQuad');
load('tNuFFT');
load('tClose');
close all
figure
plot(Nvec,tNuFFT,'ob','LineWidth',2,'HandleVisibility','off')
hold on
plot(Nvec,tNuFFT,'LineWidth',2,'DisplayName','Far-field')
plot(Nvec,tClose,'or','LineWidth',2,'HandleVisibility','off')
plot(Nvec,tClose,'LineWidth',2,'DisplayName','Close-field')
plot(Nvec,timeQuad,'oy','LineWidth',2,'HandleVisibility','off')
plot(Nvec,timeQuad,'LineWidth',2,'DisplayName','Assembling')
xlabel('Number of charges')
ylabel('Computational time (s)')
set(gca,'Fontsize',24)
box on
legend show
legend boxoff 

matlab2tikz('C:\Users\Martin\Documents\Cours\SCSD\Latex\ComptesRendus\computationalTimeLaplace.tex'); 



