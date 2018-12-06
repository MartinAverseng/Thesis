
clear all;
close all
N = 100;
A1 = -3;
B1 = 2;
A = [A1,1]; B = [B1,1];
n = [0,1];
tau = [1,0];
X = 2*4*rand(N,2)-4;
X = [X; (2*A + B)/3];
plot([A(1),B(1)],[A(2),B(2)]);
hold on
plot(X(:,1),X(:,2),'o');
[compNorm,compTang] = gradlogr(A,B,X,n,tau);

% Comparaison avec la fonction intégrale de matlab
% composante normale :
sol = compNorm; ref = 0*sol;
for i = 1:N
    Xi = X(i,:);
    d = Xi(2);
    a = Xi(1);
    fun = @(y)(-d./((y-a).^2 + d.^2));
    ref(i) = integral(fun,A1,B1);
end

norm(ref - sol,'inf')

% composante tangentielle :
sol = compTang; ref = 0*sol;
for i = 1:N
    Xi = X(i,:);
    d = Xi(2);
    a = Xi(1);
    fun = @(y)((y-a)./((y-a).^2 + d.^2));
    ref(i) = integral(fun,A1,B1);
end

norm(ref - sol,'inf')