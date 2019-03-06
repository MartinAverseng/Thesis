function [normX] = PiNormFunX(X,x,y)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin == 1
    x = 0; y = 0;
end

if size(X,2)==3
    normX = sqrt((X(:,1)-x).^2 + (X(:,2)-y).^2);
elseif size(X,1)==3
    normX = PiNormFunX(X.',x,y).';
else
    error('not possible');
end


end

