function[c] = parabola()

x = @(t)(t);
y = @(t)(t.^2 - 0.25);
I = [-1,1];
c = SimpleCurve(x,y,I);

end