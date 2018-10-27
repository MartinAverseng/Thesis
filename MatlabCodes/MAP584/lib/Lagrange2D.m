function ef = Lagrange2D(dom, ordre)
global mesh
ef.type = 'Lagrange';
ef.dim = dom.dim;
ef.label = dom.label;
ef.ordre = ordre;
ef.base = @base;
ef.ddl2int = @ddl2int;
Nvert = size(mesh.vertices,1);
Ntri = size(mesh.triangles,1);
Nedge = size(mesh.edges,1);
switch ordre
    case 0 % EF P0
       ef.Nddl = Ntri;
       ef.ddl = (mesh.vertices(mesh.triangles(:,1),:)...
           + mesh.vertices(mesh.triangles(:,2),:)...
           + mesh.vertices(mesh.triangles(:,3),:))/3;
       ef.lblDdl = mesh.lab_triangles;
    case 1 % EF P1
       ef.Nddl = Nvert;
       ef.ddl = mesh.vertices;
       ef.lblDdl = mesh.lab_vertices;
    case 2 % EF P2;
       ef.Nddl = Nvert + Nedge;
       ef.ddl = [ mesh.vertices ; ...
           (mesh.vertices(mesh.edges(:,1),:)+mesh.vertices(mesh.edges(:,2),:))/2];
       ef.lblDdl = [mesh.lab_vertices; mesh.lab_edges];
    otherwise
        error('FEspace.m : unknown order of finite element')
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function A = ddl2int(domain, idx, ef, Xhat)
global mesh
[phi, numDdl] = base(ef.ordre, Xhat);
Nddlloc = size(phi{1},2);
Nelt = size(idx,1);
Nintloc = size(Xhat,1);

% initialisation des vecteurs i,j,a
i = zeros(Nelt,Nintloc,Nddlloc); % numéro global des pts integration
j = zeros(Nelt,Nintloc,Nddlloc); % numéro global des pts ddl
for jloc = 1 : Nddlloc
    for iloc = 1 : Nintloc
        i(:,iloc,jloc) = ((0:Nelt-1))*Nintloc + iloc;
        j(:,iloc,jloc) = numDdl(idx,jloc);
    end
end

switch ef.op
    case 'Id'
        a = ones(Nelt,1)*phi{1}(:)';
    case 'Dx'
        a = (mesh.T22(idx)./mesh.detT(idx)) * phi{2}(:)' - (mesh.T21(idx)./mesh.detT(idx)) * phi{3}(:)';
    case 'Dy'
        a = -(mesh.T12(idx)./mesh.detT(idx)) * phi{2}(:)' + (mesh.T11(idx)./mesh.detT(idx)) * phi{3}(:)';
    otherwise
        error(['Non recognized operator ',op]);
end

A = sparse(i(:),j(:),a(:),Nelt*Nintloc,ef.Nddl);
end


xA = vtx(tri(:,1),1); yA = vtx(tri(:,1),2);
xB = vtx(tri(:,2),1); yB = vtx(tri(:,2),2);
xC = vtx(tri(:,3),1); yC = vtx(tri(:,3),2);

m.T11 = xB - xA; m.T12 = xC - xA;
m.T21 = yB - yA; m.T22 = yC - yA;
m.vol_triangles = m.T11.*m.T22 - m.T21.*m.T12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [phi, numDdl] = base(ordre, X)
global mesh
[N,dim]=size(X);
ONE   = ones(N,1);
ZERO  = zeros(N,1);
Nvert = size(mesh.vertices,1);
Ntri = size(mesh.triangles,1);
Nedge = size(mesh.edges,1);
if dim == 3
    x = X(:,1);
    y = X(:,2);
    z = X(:,3);
    switch ordre
        case 0 % elements P0
            phi{1} = ONE;     % values
            phi{2} = ZERO;    % x-derivatives 
            phi{3} = ZERO;    % y-derivatives
            phi{4} = ZERO;    % z-derivatives
            numDdl = (1:Ntet)';
        case 1 % elements P1
            phi{1} = [ 1-x-y-z, x   , y    , z ];
            phi{2} = [ -ONE   , ONE , ZERO , ZERO];
            phi{3} = [ -ONE   , ZERO, ONE  , ZERO];
            phi{4} = [ -ONE   , ZERO, ZERO , ONE];
            numDdl = mesh.tet;
        otherwise
            error('Lagrange2D.m : unknown Lagrange 2D order')
    end
elseif dim == 2
    x = X(:,1);
    y = X(:,2);
    switch ordre
        case 0 % elements P0
            phi{1} = ONE;     % values
            phi{2} = ZERO;    % x-derivatives 
            phi{3} = ZERO;    % y-derivatives
            numDdl = (1:Ntri)';
        case 1 % elements P1
            phi{1} = [ 1-x-y, x   , y    ];
            phi{2} = [ -ONE , ONE , ZERO ];
            phi{3} = [ -ONE , ZERO, ONE  ];
            numDdl = mesh.triangles;
        case 2 % elements P2
            phi{1} = [(1-x-y).*(1-2*x-2*y), x.*(2*x-1)   , y.*(2*y-1),  ...
                4*x.*y             , 4*y.*(1-x-y) , 4*x.*(1-x-y) ];
            phi{2} = [-3+4*(x+y)          , 4*x-1        , ZERO      ,  ...
                4*y                , -4*y         , 4*(1-2*x-y)  ];
            phi{3} = [-3+4*(x+y)          , ZERO         , 4*y-1     ,  ...
                4*x                , 4*(1-x-2*y)  , -4*x         ];
            numDdl = [mesh.triangles, Nvert + mesh.connectivity.ElementEdges];
        otherwise
            error('Lagrange2D.m : unknown Lagrange 2D order')
    end
else % dim == 1. We look for traces
    x = X(:,1);
    switch ordre
        case 0 % elements P0
            phi{1} = ONE;
            phi{2} = ZERO;
            numDdl = [];
        case 1 % elements P1
            phi{1} = [ 1-x, x];
            phi{2} = [ -ONE , ONE ];
            numDdl = mesh.edges;
        case 2 % elements P2
            phi{1}   = [(1-x).*(1-2*x), x.*(2*x-1), 4*x.*(1-x)];
            phi{2} = [4*x-3, 4*x-1, -8*x+4]; 
            numDdl = [mesh.edges, Nvert + (1:Nedge)'];
        otherwise
            error('Lagrange2D.m : unknown Lagrange 2D order')
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
