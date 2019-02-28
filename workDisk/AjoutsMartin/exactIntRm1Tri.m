function [I] = exactIntRm1Tri(A,B,C,X)
% Computes the value of int_T 1/|X - Y| dT(Y) where |.| denotes euclidean norm and 
% T is the triangle defined by the vertices A, B and C.  

S = [[A(:)'; B(:)'; C(:)'] zeros(3,1)];
X = [X(:)' 0];
M = msh(S,[1 2 3]);
n = M.nrm; % Normal vector
tau = cell2mat(M.tgt); % Tangent vectors
nu  = cell2mat(M.nrmEdg); % Normal vectors to the edges
tau  = reshape(tau,3,3)';
nu = reshape(nu,3,3)';
[I,~,~,~] = domSemiAnalyticInt3D(X,S,M.nrm,tau,nu);


end

