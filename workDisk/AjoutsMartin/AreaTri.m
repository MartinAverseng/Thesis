function [A] = AreaTri(A,B,C)


A=abs(A(:,1).*(B(:,2)-C(2))+B(:,1).*(C(:,2)-A(:,2))+C(:,1).*(A(:,2)-B(:,2)))/2.0;


end

