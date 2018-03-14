function[X,Y,V,Xaxis,Yaxis] = GridAndUniformCircleData(M,N)


Nx = round(sqrt(M));
Ny = Nx;
Xmin = -1.3;
Xmax = 1.3;

Xaxis = linspace(Xmin,Xmax,Nx);
Yaxis = linspace(Xmin,Xmax,Ny)';
gridX = repmat(Xaxis,Ny,1);
gridY = repmat(Yaxis,1,Nx);
gridX = gridX(:);
gridY = gridY(:);

X = [gridX,gridY];
r = 0.5;
perturbation = rand(N,1)*0;
theta = (linspace(0,1,N)'+perturbation)*2*pi;
Y(:,1) = r*cos(theta);
Y(:,2) = r*sin(theta);


phases = (randn(N,1)+perturbation)*2*pi;
V = (1+perturbation).*exp(1i*phases);
V = V/norm(V,1);


end