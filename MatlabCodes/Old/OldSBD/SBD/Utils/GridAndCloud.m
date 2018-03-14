function[X,Y,V,Xaxis,Yaxis] = GridAndCloud(Ngrid,Ncharges)

Nx = round(sqrt(Ngrid));
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
Y = 2*rand(Ncharges,2)-1;
V = 2*rand(Ncharges,1)-1;
V = V/norm(V,1);

end