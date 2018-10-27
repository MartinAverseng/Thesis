function [ c , incWave] = unitSegment(k)

c = openline(-1,1);



if and(nargout == 2,nargin ==1)
    theta_inc = pi/2;
    X = R2toRfunc.X; Y = R2toRfunc.Y;
    incWave = exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));
end


end
