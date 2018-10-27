function ef = Lagrange(dom, ordre)
global mesh
ef.type = 'Lagrange';
ef.dim = dom.dim;
ef.label = dom.label;
ef.ordre = ordre;
ef.base = @base;
ef.ddl2int = @ddl2int;

Nvert = size(mesh.vertices,1);
Nedge = size(mesh.edges,1);
Ntri = size(mesh.triangles,1);
Ntet = size(mesh.tetras,1);
switch dom.dim
    case 1
        idx = find(ismember(mesh.edges,dom.label));
        switch ordre
            case 0 % EF P0
                ef.Nddl = size(idx,1);
                ef.ddl = (mesh.vertices(mesh.edges(idx,1),:)...
                    + mesh.vertices(mesh.edges(idx,2),:))/2;
                ef.lblDdl = mesh.lab_edges(idx);
            case 1 % EF P1
                % number of vertices in the domain
                elts = mesh.edges(idx,:);
                vertices = unique(elts(:));
                Nvert = size(vertices);
                ef.Nddl = Nvert;
                ef.ddl = mesh.vertices(vertices(:),:);
                ef.lblDdl = mesh.lab_vertices(vertices(:));
            case 2 % EF P2;
                % number of vertices in the domain
                elts = mesh.edges(idx,:);
                vertices = unique(elts(:));
                Nvert = size(vertices);
                Nedge = size(idx,1);
                ef.Nddl = Nvert + Nedge;
                ef.ddl = [ mesh.vertices(vertices,:) ; ...
                    (mesh.vertices(mesh.edges(idx,1),:)+mesh.vertices(mesh.edges(idx,2),:))/2];
                ef.lblDdl = [mesh.lab_vertices(vertices(:)); mesh.lab_edges(idx)];
            otherwise
                error('FEspace.m : unknown order of finite element')
        end        
    case 2
        idx = find(ismember(mesh.triangles,dom.label));
        switch ordre
            case 0 % EF P0
                ef.Nddl = size(idx,1);
                ef.ddl = (mesh.vertices(mesh.triangles(idx,1),:)...
                    + mesh.vertices(mesh.triangles(idx,2),:)...
                    + mesh.vertices(mesh.triangles(idx,3),:))/3;
                ef.lblDdl = mesh.lab_triangles(idx);
            case 1 % EF P1
                % number of vertices in the domain
                elts = mesh.triangles(idx,:);
                vertices = unique(elts(:));
                Nvert = size(vertices);
                ef.Nddl = Nvert;
                ef.ddl = mesh.vertices(vertices,:);
                ef.lblDdl = mesh.lab_vertices(vertices);
            case 2 % EF P2;
                % number of vertices in the domain
                elts = mesh.triangles(idx,:);
                vertices = unique(elts(:));
                Nvert = size(vertices);

                ef.Nddl = Nvert + Nedge;
                ef.ddl = [ mesh.vertices(vertices,:) ; ...
                    (mesh.vertices(mesh.edges(:,1),:)+mesh.vertices(mesh.edges(:,2),:))/2];
                ef.lblDdl = [mesh.lab_vertices; mesh.lab_edges];
            otherwise
                error('FEspace.m : unknown order of finite element')
        end
    case 3
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function A = ddl2int(domain, idx, Xhat, ef)
global mesh
phi = base(ef.ordre, Xhat);
Nddlloc = size(phi{1},2);
Nelt = size(idx,1);
Nintloc = size(Xhat,1);

%%%%% Il faut calculer numDdl et le nombre total de Ddl Nddl!!!!!!!!!!

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
function phi = base(ordre, X)
% Computes the basis functions in the reference element
global mesh
[N,dim]=size(X);
ONE   = ones(N,1);
ZERO  = zeros(N,1);
switch dim
    case 3
        x = X(:,1);
        y = X(:,2);
        z = X(:,3);
        switch ordre
            case 0 % elements P0
                phi{1} = ONE;     % values
                phi{2} = ZERO;    % x-derivatives
                phi{3} = ZERO;    % y-derivatives
                phi{4} = ZERO;    % z-derivatives
            case 1 % elements P1
                phi{1} = [ 1-x-y-z, x   , y    , z ];
                phi{2} = [ -ONE   , ONE , ZERO , ZERO];
                phi{3} = [ -ONE   , ZERO, ONE  , ZERO];
                phi{4} = [ -ONE   , ZERO, ZERO , ONE];
            otherwise
                error('Lagrange.m : unknown Lagrange 3D order')
        end
    case 2
        x = X(:,1);
        y = X(:,2);
        switch ordre
            case 0 % elements P0
                phi{1} = ONE;     % values
                phi{2} = ZERO;    % x-derivatives
                phi{3} = ZERO;    % y-derivatives
            case 1 % elements P1
                phi{1} = [ 1-x-y, x   , y    ];
                phi{2} = [ -ONE , ONE , ZERO ];
                phi{3} = [ -ONE , ZERO, ONE  ];
            case 2 % elements P2
                phi{1} = [(1-x-y).*(1-2*x-2*y), x.*(2*x-1)   , y.*(2*y-1),  ...
                    4*x.*y             , 4*y.*(1-x-y) , 4*x.*(1-x-y) ];
                phi{2} = [-3+4*(x+y)          , 4*x-1        , ZERO      ,  ...
                    4*y                , -4*y         , 4*(1-2*x-y)  ];
                phi{3} = [-3+4*(x+y)          , ZERO         , 4*y-1     ,  ...
                    4*x                , 4*(1-x-2*y)  , -4*x         ];
            otherwise
                error('Lagrange.m : unknown Lagrange 2D order')
        end
    case 1
        x = X(:,1);
        switch ordre
            case 0 % elements P0
                phi{1} = ONE;
                phi{2} = ZERO;
            case 1 % elements P1
                phi{1} = [ 1-x, x];
                phi{2} = [ -ONE , ONE ];
            case 2 % elements P2
                phi{1}   = [(1-x).*(1-2*x), x.*(2*x-1), 4*x.*(1-x)];
                phi{2} = [4*x-3, 4*x-1, -8*x+4];
            otherwise
                error('Lagrange.m : unknown Lagrange 1D order')
        end
end
end
