% Validation of the mass matrix and single layer potential with eigenvalues

N   = 10;
tol = 1e-4;
typ = 'P1';

% Mesh the disk and the half sphere S2
Disku = mshMyDisk(N,1);

% Ad hoc version

S2 = Disku;
r = sqrt(sum(S2.vtx.^2,2));
theta = pi/2*(1-r);
S2.vtx(2:end,1) = S2.vtx(2:end,1)./r(2:end).*cos(theta(2:end));
S2.vtx(2:end,2) = S2.vtx(2:end,2)./r(2:end).*cos(theta(2:end));
S2.vtx(:,3) = sin(theta(:));
Disk = S2;
Disk.vtx(:,3) = 0;
Disk2 = Disk;
Disk2.wgt = S2.ndv./Disk.ndv;

% Check weighted integral 
surface = sum(Disk2.ndv);
% Domain
sigmaw = weightedDom(Disk2, 2);

% Finite elements on the sphere
Vh = fem(Disk2,typ);

% Incident wave
PW = @(X) ones(size(X,1),1);

% Mass matrix
tic
Mw = integral(sigmaw,Vh,Vh);
% Mw = (Mw + Mw')/2;
tMw = toc;

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' G(x,y) psi(y) dx dy 
tic
Sw = 1/(4*pi) .* integral(sigmaw,sigmaw,Vh,Gxy,Vh);
Srw  = 1/(4*pi) .* projRegularize(sigmaw, sigmaw, Vh, '[1/r]', Vh);
Sw = Sw + Srw;
tSw = toc

[P,D] = eig(Mw\Sw);
d = sort(diag(D),'descend');

lambdaTheo = [];
for l = 0:10
    for m = -l:l
        if mod(l + m,2) == 0
            lpm = (l + m)/2;
            lmm = (l - m)/2;
            lambdaTheo(end+1,1) = 1/4*gamma(lpm + 1/2)*gamma(lmm + 1/2)/(gamma(lpm + 1)*gamma(lmm +1));
        end
    end
end
lambdaTheo = sort(lambdaTheo,'descend');
disp([lambdaTheo(1:20) d(1:20) lambdaTheo(1:20)-d(1:20)])
% Weighted Delta matrix 
weight = @(X) (1.0000000001 - X(:,1).^2 - X(:,2).^2);
tic
Delta = integral(sigmaw,grad(Vh),weight,grad(Vh));
Delta = (Delta + Delta')/2;
tDelta = toc

[P,D] = eig(full(Mw\Delta));
d = sort(diag(D));

lambdaTheo = [];
for l = 0:10
    for m = -l:l
        if mod(l + m,2) == 0
            lpm = (l + m)/2;
            lmm = (l - m)/2;
            lambdaTheo(end+1,1) = l*(l+1) - m^2;
        end
    end
end
lambdaTheo = sort(lambdaTheo);
disp([lambdaTheo(1:20) d(1:20) lambdaTheo(1:20)-d(1:20)])


[P,lambdasProd] = eig(full(Mw^(-1)*Sw*Mw^(-1)*(Delta + 5*Mw)*Mw^(-1)*Sw));


disp('Problème ?');

