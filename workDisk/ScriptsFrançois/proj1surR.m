function Gxy = proj1surR(X,Y)

% Projected distances
Rxy = sqrt((X(:,1) - Y(:,1)).^2 + (X(:,2) - Y(:,2)).^2 );
Gxy = 1./Rxy;   

% Singularity
Gxy(Rxy<1e-12) = 0;
end

