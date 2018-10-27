T = 100;
t = linspace(0,1,T);
N = 5;
wt = zeros(T,1);
for j = 1:T
    wt(j) = wronkianBessel(N,t(j));
end
plot(t,wt);