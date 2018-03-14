%% Circular quadrature for the bessel function :

% It is trivial, the quadrature must be chosen periodic.
clear all;
close all;
x = 1000.11; y = 1000;

Ns = fix(sqrt(x^2+y^2)*0.5):fix(sqrt(x^2+y^2)*10);


for i = 1:length(Ns)
    N = Ns(i);
    ksi = [cos((0:N-1)'*2*pi/N) sin((0:N-1)'*2*pi/N)];
    quadC = mean(exp(1i*ksi*[x;y]));
    err(i) = log(abs(besselj(0,sqrt(x^2 + y^2))-quadC));
end

plot(err);
