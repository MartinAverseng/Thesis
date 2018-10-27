function[c] = ellipse()

a = 1;
b = 1;

x = @(t)(a*cos(1.4*pi/2*t));
y = @(t)(b*sin(1.4*pi/2*t));
I = [-1 1];

c = SimpleCurve(x,y,I);

end