function [ ] = surfScatter( x,y,z )


%% Little triangles
% The solution is to use Delaunay triangulation. Let's look at some
% info about the "tri" variable.

tri = delaunay(x,y);

%% Plot it with TRISURF

trisurf(tri, x, y, z);
axis vis3d

%% Clean it up

axis off
light('Position',[-50 -15 29]);
lighting phong
shading interp
colorbar EastOutside


end

