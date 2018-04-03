function [arc ] = spirale()

x = @(s)(exp(s).*cos(10*s));
y = @(s)(exp(s).*sin(10*s));
I = [-1,1];
arc = SimpleCurve(x,y,I);


end

