function plotMesh(varargin)
%
% plotmesh(varargin)
%
% varargin = m 
% => trace le maillage m
%
% varagin = "m , option"
% => trace le maillage m
% si option = 'v' ou 'V', indique respectivement les numéros des sommets ou leur label
%           = 't' ou 'T', indique respectivement les numéros des triangles ou leur label
%           = 'e' ou 'E', indique respectivement les numéros des arêtes ou leur label


if (nargin==0 || nargin>2 )
    error('plotMesh.m : Use plotMesh(m) ou plotMesh(m,option)')
end

% tracé du maillage
m=varargin{1};

vtx=m.vertices;
vtxlab=m.lab_vertices;
tri=m.triangles;
trilab=m.lab_triangles;
edg=m.edges;
edglab=m.lab_edges;
x=vtx(:,1);
y=vtx(:,2);
nbNodes=size(vtx,1);
nbTri=size(tri,1);
nbEdge=size(edg,1);


trimesh(tri,x,y,'Color','b');

if nargin == 2
    option=varargin{2};
    if findstr(option,'v')
        txt=num2str([1:nbNodes]');
        text(x,y,txt);
    end
    if findstr(option,'t')
        txt=num2str([1:nbTri]');
        xtri=(vtx(tri(:,1),1)+vtx(tri(:,2),1)+vtx(tri(:,3),1))/3;
        ytri=(vtx(tri(:,1),2)+vtx(tri(:,2),2)+vtx(tri(:,3),2))/3;
        text(xtri,ytri,txt);
    end
    if findstr(option,'e')
        txt=num2str([1:nbEdge]');
        xedg=(vtx(edg(:,1),1)+vtx(edg(:,2),1))/2;
        yedg=(vtx(edg(:,1),2)+vtx(edg(:,2),2))/2;
        text(xedg,yedg,txt);
    end
    if findstr(option,'V')
        txt=num2str(vtxlab);
        text(x,y,txt);
    end
    if findstr(option,'T')
        txt=num2str(trilab);
        xtri=(vtx(tri(:,1),1)+vtx(tri(:,2),1)+vtx(tri(:,3),1))/3;
        ytri=(vtx(tri(:,1),2)+vtx(tri(:,2),2)+vtx(tri(:,3),2))/3;
        text(xtri,ytri,txt);
    end
    if findstr(option,'E')
        txt=num2str(edglab);
        xedg=(vtx(edg(:,1),1)+vtx(edg(:,2),1))/2;
        yedg=(vtx(edg(:,1),2)+vtx(edg(:,2),2))/2;
        text(xedg,yedg,txt);
    end
end

axis equal