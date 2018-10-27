function [ L ] = secMembre(f,X)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% N = length(X);
% T = acos(X);
% T = T(end:-1:1);
% Delta = diff(T);
% phi = T(1:end-1) + Delta/2;
% 
% % Le vecteur L est alors tout simplement
% 
% L = Delta.*f(cos(phi));
% L = L(:);

N = length(X);
T = acos(X);
T = T(end:-1:1);
[phi,w] = gauss_lobatto(-1,1);
M = length(phi);
Phi = repmat(phi,1,N-1);
alpha = diff(T)/2;
beta = alpha + T(2:end);
beta = beta(:)';
Phi_k = Phi*diag(alpha)+repmat(beta,M,1);
W = repmat(w(:),1,N-1)*diag(alpha);

L = sum(f(cos(Phi_k)).*W,1)';

end

