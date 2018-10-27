function m=mesh2D(g,aireMax)

m=triangle(g,aireMax);
m.tetra=[];
m.lab_tetra=[];
m.vol_tetra = [];

xA = m.vertices(m.edges(:,1),:); xB = m.vertices(m.edges(:,2),:);
m.vol_edges = sqrt(sum((xB-xA).^2,2));

vtx = m.vertices;
tri = m.triangles;
xA = vtx(tri(:,1),1); yA = vtx(tri(:,1),2);
xB = vtx(tri(:,2),1); yB = vtx(tri(:,2),2);
xC = vtx(tri(:,3),1); yC = vtx(tri(:,3),2);

m.T11 = xB - xA; m.T12 = xC - xA;
m.T21 = yB - yA; m.T22 = yC - yA;
m.vol_triangles = m.T11.*m.T22 - m.T21.*m.T12;
