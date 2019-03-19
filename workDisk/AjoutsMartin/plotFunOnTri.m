function[] = plotFunOnTri(A,B,C,f)

N1 = @(x,y)(1 - x - y);
N2 = @(x,y)(x);
N3 = @(x,y)(y);

P = @(x,y)(A(1)*N1(x,y) + B(1)*N2(x,y) + C(1)*N3(x,y));
Q = @(x,y)(A(2)*N1(x,y) + B(2)*N2(x,y) + C(2)*N3(x,y));

% Points on the reference triangle:
[x,y] = meshgrid(0:0.01:1);
x(x+y>1) = nan;
y(x+y>1) = nan;

z = f(P(x,y),Q(x,y));
z(x + y>1) = nan;
surf(P(x,y),Q(x,y),z);

end