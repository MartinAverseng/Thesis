% Tests3
% On va tester les formules de quadratures sur la sphère. 

N = 24;
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
figure
plot(Disk2);

Disk2.wgt = S2.ndv./Disk.ndv;
% Nrm = S2.nrm;
%Disk2.alpha = sum(S2.vtx(S2.elt(:,1),:).*Nrm,2);
% Check weighted integral 
surface = sum(Disk2.ndv);
disp(abs(surface - 2*pi));

hold on;


I = 0;
for i = 1:size(Disk2.elt,1)
    tri = Disk2.elt(i,:);
    triCoords = Disk2.vtx(tri,:);
    A = triCoords(1,:);
    B = triCoords(2,:);
    C = triCoords(3,:);
    [X,W] = quadSphTri(A,B,C,2);
    I = I + sum(W);
end
disp(abs(I - 2*pi))

sigma = weightedDom(Disk2,2);
integral(
