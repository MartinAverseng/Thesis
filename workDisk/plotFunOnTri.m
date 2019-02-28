function[] = plotFunOnTri(A,B,C,f)

N1 = @(x,y)(1 - x - y);
N2 = @(x,y)(x);
N3 = @(x,y)(y);

P = @(x,y)(A(1)*N1(x,y) + B(1)*N2(x,y) + C(1)*N3(x,y));
Q = @(x,y)(A(2)*N1(x,y) + B(2)*N2(x,y) + C(2)*N3(x,y));

% Points on the reference triangle:
x = linspace(0,1,100);
y = linspace(0,1,100);
X = [];
Y = [];
for i = 1:length(x)
    yi = y(y+x(i)<=1);
    xi = x(i)*ones(1,length(yi));
    X = [X xi];
    Y = [Y yi];
end

x = P(X,Y);
y = Q(X,Y);

scatter3(x,y,f(x',y'));

end