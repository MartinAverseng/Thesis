function [I] = gaussQuadTri(freg,A,B,C,N)
% This function evaluates \iint_K f(x,y) dxdy using
% the Gaussian quadrature of order N where K is a
% triangle with vertices (x1,y1), (x2,y2) and (x3,y3).
x1 = A(1); y1 = A(2);
x2 = B(1); y2 = B(2);
x3 = C(1); y3 = C(2);
xw = TriGaussPoints(N); % get quadrature points and weights
% calculate the area of the triangle
A=abs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/2.0;
% find number of Gauss points
NP=length(xw(:,1));
z = 0.0;
for j = 1:NP
x = x1*(1-xw(j,1)-xw(j,2))+x2*xw(j,1)+x3*xw(j,2);
y = y1*(1-xw(j,1)-xw(j,2))+y2*xw(j,1)+y3*xw(j,2);
z = z + freg(x,y)*xw(j,3);
end
I = A*z;

return
end
