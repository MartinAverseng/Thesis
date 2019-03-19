function [I] = explicitIntRegularize(A,B,C,theta)
% Attemps to compute int_{Tri} f(s,t) ds dt
% where f = @(s,t)(1./sqrt(t.^2 + s.^2*cos(theta)^2))
% the triangle Tri is defined by its vertex A,B,C in R^2

% Matlab answer :
f = @(s,t)(1./sqrt(t.^2 + s.^2*cos(theta)^2));
I = integral2tri(f,A,B,C);


end

