%% Validation of bessel zeros search 

close all;

%% Dirichlet condition

N = 30; % Number of zeros we search


subplot(1,2,1)
D = 2; % Dimension 2
k = D/2-1;
zs = BesselZeros(Inf,k,N);
t = linspace(0,zs(end)+1,1000);
plot(t,besselj(k,t));
hold on
plot(zs,zeros(size(zs)),'.','MarkerSize',15)
title(sprintf('Dirichlet D=%d',D))

subplot(1,2,2)
D = 3; % Dimension 3
k = D/2-1;
zs = BesselZeros(Inf,k,N);
t = linspace(0,zs(end)+1,1000);
plot(t,besselj(k,t));
hold on
plot(zs,zeros(size(zs)),'.','MarkerSize',15)
title(sprintf('Dirichlet D=%d',D))

ns = zs/pi;
assert(norm(ns-(0:N-1))<1e-5);


%% Robin condition 


N = 50; % Number of zeros we search
c = 1;
figure
subplot(1,2,1)
D = 2; % Dimension 2
k = D/2-1;
zs = BesselZeros(c,k,N);
t = linspace(0,zs(end)+1,1000);
plot(t,c*besselj(k,t)+0.5*t.*(besselj(k-1,t)-besselj(k+1,t)));
hold on
plot(zs,zeros(size(zs)),'.','MarkerSize',15)
title(sprintf('Robin D=%d',D))

subplot(1,2,2)
D = 3; % Dimension 3
k = D/2-1;
zs = BesselZeros(c,k,N);
t = linspace(0,zs(end)+1,1000);
plot(t,c*besselj(k,t)+0.5*t.*(besselj(k-1,t)-besselj(k+1,t)));
hold on
plot(zs,zeros(size(zs)),'.','MarkerSize',15)
title(sprintf('Robin D=%d',D))

