function [ rho ] = BesselZeros(N,startFreq)

f  = @(r) besselj(0,r);
df = @(r) - besselj(1,r);
init_guess = startFreq + (0:N-1)'*pi;
rho = newton(f,df,init_guess);


end

