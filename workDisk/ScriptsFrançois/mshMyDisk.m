function mesh = mshMyDisk(N,rad)

% Radial discretisation
dr = rad/N;
r  = dr:dr:rad;

% Angular uniform discretization
rho = cell(length(r),1); theta = rho;
for ir = 1:length(r)
    dtheta = 2*pi/(6*ir);
    theta{ir} = (0:dtheta:2*pi-dtheta)';
    rho{ir}   = r(ir)*ones(length(theta{ir}),1);
end

% Carthesian coordinates
[x,y] = pol2cart(cell2mat(theta),cell2mat(rho));
x = [0;x];
y = [0;y];
   
% Delaunay triangulation
DT = delaunayTriangulation(x,y);

% Final mesh
elt  = DT.ConnectivityList;
vtx  = [DT.Points,zeros(size(DT.Points,1),1)];
mesh = msh(vtx,elt);
end
