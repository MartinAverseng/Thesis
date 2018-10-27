%% Décomposition de l'unité en sinus cardinal :

clear all;
close all;

eps = 1e-3; % target error threshold (L2)
rho = 0.3; % Range of close interaction
P = fix(SCSD_orderEst(rho,eps))+1;

fprintf('Model order used : %d\n',P)

[Beta,err1] = SCSDcoeff(rho,P);


for k = 1:3
    [q{k}] = LanczosCoeff(k,P);
    leg{k} = sprintf('k=%d',k);
end

q{k+1} = Beta;
leg{k+1} = 'SCSD';

Npoint = 100;
[h1,h2] = compareQuads(q,leg,Npoint,rho,eps);

saveas(h1,'Quadratures/figures/CompareQuadNormal');
saveas(h2,'Quadratures/figures/CompareQuadLog');

rhos = [0.5; 0.1; 0.01];
kmax = 3;
eps = 1e-6;
PP = zeros(kmax+2,length(rhos));
for i = 1:length(rhos)
    PP(1:kmax,i) = reachAccuracyLanczos(rhos(i),kmax,Npoint,eps);
end
for i = 1:length(rhos)
    PP(kmax+1,i) = reachAccuracySCSD(rhos(i),Npoint,eps);
    PP(kmax+2,i) = SCSD_orderEst(rhos(i),eps);
end
