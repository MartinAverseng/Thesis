close all;
addpath('lib');
g.points = [0 0; 2 0;2 1; 1 1; 1 2; 0 2;];
g.segments =  [1 2; 2 3; 3 4; 4 5; 5 6; 6 1];
h = 0.1;
m = triangle(g,0.023);
plotMesh(m);
ver = m.vertices;
tri = m.triangles;
N = size(tri,1); % Nombre de triangle
areasRefine = zeros(N,1);
for k = 1:N
    trik = tri(k,:)';
    versk = ver(trik,:);
    gravk = mean(versk);
    areasRefine(k) = h*norm(gravk - [1,1])^2;
end

m2 = refine(m,areasRefine);