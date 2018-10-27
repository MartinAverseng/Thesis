function [arc,incWave ] = spirale(k)

a= 1;
b= 0.2;
c= 10;
s0 = -0.7;
x = @(s)(a*exp(b*c*(s+s0)).*cos(c*(s+s0)));
y = @(s)(a*exp(b*c*(s+s0)).*sin(c*(s+s0)));
I = [-1,1];
arc = SimpleCurve(x,y,I);

if and(nargout == 2,nargin ==1)
    theta_inc = 0 + pi/2;
    X = R2toRfunc.X; Y = R2toRfunc.Y;
    incWave = exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));
end



end

