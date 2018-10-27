function [ft] = f_n(t,H)

H = H(:);
N = (length(H)-1)/2;
vec = (-N:N)';
ft = (H'*exp(1i*vec*t)).*exp(1i*cos(t));



end

