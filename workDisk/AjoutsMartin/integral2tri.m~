function [I] = integral2tri(f,A,B,C)
% Computes the approximation by adaptive integration of int(f) on ABC. 
AT = abs(A(1)*(B(2) - C(2))  + B(1)*(C(2) - A(2)) + C(1)*(A(2) - B(2)))/2;

N1 = @(x,y)(1 - x - y);
N2 = @(x,y)(x);
N3 = @(x,y)(y);

P = @(x,y)(A(1)*N1(x,y) + B(1)*N2(x,y) + C(1)*N3(x,y));
Q = @(x,y)(A(2)*N1(x,y) + B(2)*N2(x,y) + C(2)*N3(x,y));

xmin = 0; xmax = 1;
ymin = 0; ymax = @(x)(1 - x);

F = @(x,y)(f(P(x,y),Q(x,y)));

I = 2*AT*integral2(F,xmin,xmax,ymin,ymax);



end

