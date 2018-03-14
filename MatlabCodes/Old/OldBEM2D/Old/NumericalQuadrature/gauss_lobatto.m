function [x,w] = gauss_lobatto(a,b)
% Rescale les poids du segment [-1,1] sur [a,b]

if nargin ==0
    a = -1;
    b = 1;
end

alpha = (b-a)/2;
beta = (a+b)/2;

a = sqrt(1/21*(7-2*sqrt(7)));
wa = 1/30*(14+sqrt(7));
wb = 1/30*(14-sqrt(7));
w1 = 1/15;
b = sqrt(1/21*(7+2*sqrt(7)));
x0 = [-1;-b;-a;a;b;1];
w0 = [w1;wb;wa;wa;wb;w1];

x = alpha*x0 + beta;
w = alpha*w0;



end

