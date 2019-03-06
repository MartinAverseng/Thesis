function [xw] = TriGaussABC(A,B,C,N)
% Computes the Gaussian quadrature of order N where A B C is a triange in
% R2.

x1 = A(1); y1 = A(2);
x2 = B(1); y2 = B(2);
x3 = C(1); y3 = C(2);
xwref = TriGaussPoints(N); % get quadrature points and weights
xw = 0*xwref;
% calculate the area of the triangle
A=abs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/2.0;
% find number of Gauss points

xw(:,1) = x1*(1-xwref(:,1)-xwref(:,2))+x2*xwref(:,1)+x3*xwref(:,2);
xw(:,2) = y1*(1-xwref(:,1)-xwref(:,2))+y2*xwref(:,1)+y3*xwref(:,2);
xw(:,3) = A*xwref(:,3);

end