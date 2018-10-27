%
% triangle : 2D triangle mesh generator
%
%   mesh = triangle(geom, aire);
%
%   Mesh a given 2D domain in triangles
%
% Inputs :
%         geom   structure describing the geometry (see below)
%         aire   maximal area of a triangle 
%
% Output :
%         mesh   structure describing the generated mesh
%
% Geometry structure :
%
%         points        array nP x 2 of coordinates of the vertices of the boundary
%                         (either interior or exterior)
%         segments      array nS x 2 of boundary segments (either interior or exterior)
%         lab_points    (optional) array nP x 1 of points labels
%         lab_segments  (optional) array nS x 1 of segments labels
%         holes         (optional) array nH x 2 of points inside each hole
%         regions       (optional) array nR x 2 of points inside each region
%         lab_regions   (optional) array nR x 1 of region's labels
%
% Mesh structure :
%
%         vertices      array nV x 2 of vertices coordinates
%         triangles     array nT x 3 of vertices numers associated to each triangle
%         edges         array nE x 2 of vertices associated to each edge
%         lab_vertices  array nV x 1 of vertices labels
%         lab_triangles array nT x 1 of triangles labels
%         lab_edges     array nE x 1 of edge labels
%         
%         connectivities structure of mesh connectivity
%            ElementNodes : for each triangle, its vertices
%            ElementEdges : for each triangle, its edges
%            NodeElements : for each vertex, a list of triangles containing this vertex
%            NodeEdges    : for each vertex, the edges based on that vertex
%            EdgeElements : for each edge, the 2 neighboring elements
%            EdgeNodes    : for each edge the 2 vertices
%
%     Example :
%         g.points = [ 0. 0. ; 1. 0. ; ... ];
%         g.segments = [ 1 2 ; 2 3 ; 3 4 ; ... ];
%         m = triangle(g, 0.01)
%
%    (c)   Marc Tajchman (2002)
