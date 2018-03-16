function [arc ] = spirale()

x = @(s)(exp(s).*cos(5*s));
y = @(s)(exp(s).*sin(5*s));
I = [-1,1];
arc = SimpleCurve(x,y,I);


end

