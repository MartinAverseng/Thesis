function [nodes] = tchebNodes( N )
% Points de Tchebitchev définis par
angles = anglesTcheb(N);
nodes = cos(angles);

% Ces points sont le projeté sur (-1,1) d'un ensemble de points bien
% répartis sur le cercle unité de R².


end

