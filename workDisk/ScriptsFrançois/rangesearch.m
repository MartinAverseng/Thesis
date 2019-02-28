function [IX,IY,A] = myRangeSearch(X,Y,rho)
%
% Compute the interactions between the sets of points X and Y that are
% smaller than a given bound rho.
% Three calls are available:
% - C = myRangeSearch(X,Y,rho)
%   Here the returned sparse matrix C contains the distances.
%   In other words C_ij = |X_i-Y_j|
% - [I,D] = myRangeSearch(X,Y,rho)
%   This is the format compatible with Matlab's routine rangesearch 
%   (in the statistical toolbox).
% - [IX,IY,D] = myRangeSearch(X,Y,rho)
%   Here, three vectors are returned. I that contains the indices
%   in the table X, IY those in Y and D the corresponding distances.
%
% (c) François Alouges 2018
%

% Build the englobing box
NX = size(X,1);
NY = size(Y,1);
dim = size(X,2);
M = max([X;Y],[],1)+rho/10;
m = min([X;Y],[],1)-rho/10;

% Divide the box into cells of size rho
L = ceil((M-m)/rho);
Nbox = prod(L);
CPL = cumprod([1 L])';

% Locate the box for each point in X and build the correspondance matrix
locX = floor((X - ones(NX,1)*m)./rho);
numX = ones(NX,1) + locX*CPL(1:dim);
Aleft = sparse(numX, 1:NX, ones(NX,1), Nbox, NX);

% Same for Y
locY = floor((Y - ones(NY,1)*m)./rho);
numY = ones(NY,1) + locY*CPL(1:dim);
Aright = sparse(numY, 1:NY, ones(NY,1), Nbox, NY);

% Compute the shifts for the bounding boxes
diags = 0;
for i=1:dim
    diags = [diags-CPL(i) diags diags+CPL(i)];
end
Ndiags = size(diags,2);

% Build the inner sparse matrix
boxY = unique(numY);
boxX = unique(numX);
idxtot = cell(Ndiags,1);
jdxtot = cell(Ndiags,1);
coeftot  = cell(Ndiags,1);
idxtot = cell(Ndiags,1);
jdxtot = cell(Ndiags,1);
coeftot = cell(Ndiags,1);

for i=1:Ndiags
    jdx = boxY;
    idx = jdx + diags(i);
    test = find(ismember(idx,boxX));
    B = sparse(idx(test),jdx(test),ones(size(test)),Nbox,Nbox);

    % Build the interaction matrix
    C = Aleft'*B*Aright;

    % Erase the distances that are too big and prepare output
    [idx,jdx,~] = find(C);
%    dist = sum((X(idx(:),:)-Y(jdx(:),:)).^2,2);
    dist = (X(idx(:),1)-Y(jdx(:),1)).^2 + (X(idx(:),2)-Y(jdx(:),2)).^2 + (X(idx(:),3)-Y(jdx(:),3)).^2;
    test = (dist <= rho^2);
    jdxtot{i} = jdx(test);
    idxtot{i} = idx(test);
    coeftot{i} = sqrt(dist(test));
end
jdxtot = cell2mat(jdxtot);
idxtot = cell2mat(idxtot);
coeftot = cell2mat(coeftot);

[jdx,I]=sort(jdxtot);
idx = idxtot(I);
coef = coeftot(I);

% Output as array of cells like Matlab's rangesearch
if nargout == 3
    IX = idx;
    IY = jdx;
    A = coef;
elseif nargout == 2
    % Matlab Rangesearch Format    
    nbr = accumarray(jdx,1,[NY,1]);  
    % Final output
    IX = cell(NY,1);
    IY = cell(NY,1);
    istart = 1;
    for n = 1:NY
        if nbr(n) == 0
            IX{n} = zeros(1,0);
            IY{n} = zeros(1,0);
        else
            IX{n} = idx(istart:istart+nbr(n)-1)';
            IY{n} = coef(istart:istart+nbr(n)-1)';
        end
        istart = istart + nbr(n);
    end
else
    IX = sparse(idx, jdx, coef, NX, NY);
end
end