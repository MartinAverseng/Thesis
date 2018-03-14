function[X,Y,V] = randomCloud(M,N)
X = randn(M,2);
Y = randn(N,2);
V = randn(N,1);
V = V/norm(V,1);
end