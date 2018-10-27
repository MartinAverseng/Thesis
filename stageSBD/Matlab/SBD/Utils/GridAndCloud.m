function[X,Y,V,Xaxis,Yaxis] = GridAndCloud(M,N)

Nx = round(sqrt(M));
Ny = Nx;
Xmin = -1;
Xmax = 1;

Xaxis = linspace(Xmin,Xmax,Nx);
Yaxis = linspace(Xmin,Xmax,Ny)';
gridX = repmat(Xaxis,Ny,1);
gridY = repmat(Yaxis,1,Nx);
gridX = gridX(:);
gridY = gridY(:);

X = [gridX,gridY];
Y = randn(N,2);
V = randn(N,1);
V = V/norm(V,1);

end