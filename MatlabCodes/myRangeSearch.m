function C = myRangeSearch(rho,varargin)
%
% Appel C = myRangeSearch(rho,X) ou C = myRangeSearch(rho,X,Y)
%
% Construit la matrice creuse C_ij =|X_i-X_j] ou [X_i-Y_j]
% pour les couples de points dont la distance est plus petite que rho.
%
% (c) François Alouges 2018
%
X = varargin{1};
NX = size(X,1);
dim = size(X,2);
M = max(X,[],1)+rho/2;
m = min(X,[],1)-rho/2;
if nargin==3
    Y = varargin{2};
    NY = size(Y,1);
    MY = max(Y,[],1)+rho/2;
    mY = min(Y,[],1)-rho/2;
    M = max(M,MY);
    m = min(m,mY);
else
    Y = X;
    NY = NX;
end
L = ceil((M-m)/rho);
NN = prod(L);
CPL = cumprod([1 L])';
locX = floor((X - ones(NX,1)*m)./rho);
num = ones(NX,1) + locX*CPL(1:dim);
Aleft = sparse(num, 1:NX, ones(NX,1), NN, NX);
if nargin==3
    locY = floor((Y - ones(NY,1)*m)./rho);
    num = ones(NY,1) + locY*CPL(1:dim);
    Aright = sparse(num, 1:NY, ones(NY,1), NN, NY);
else
    Aright = Aleft;
end
% Matrice 3, 9 ou 27 diagonales
if dim==1
    diags = 1;
elseif dim==2
    diags = [1 L(1)-1 L(1) L(1)+1];
elseif dim==3
    diags = [1 L(1)-1 L(1) L(1)+1  ...
        L(1)*(L(2)-1)-1 L(1)*(L(2)-1) L(1)*(L(2)-1)+1 ...
        L(1)*L(2)-1 L(1)*L(2) L(1)*L(2)+1 ...
        L(1)*(L(2)+1)-1 L(1)*(L(2)+1) L(1)*(L(2)+1)+1];
else
    error('dimension not supported')
end
B = spdiags(ones(NN,27),sort([-diags 0 diags]),NN,NN);
C = Aleft'*B*Aright;
[i,j,~] = find(C);
dist = sum((X(i(:),:)-Y(j(:),:)).^2,2);
idx = find(dist <= rho^2);
C = sparse(i(idx),j(idx),sqrt(dist(idx)),NX,NY);
end