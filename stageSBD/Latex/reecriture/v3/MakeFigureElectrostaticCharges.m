clear all
close all

Ncharges = 10000;
Ngrid = 1000000;
[X,Y,V,Xaxis,Yaxis] = GridAndCloud(Ngrid,Ncharges);
X = X/(2*sqrt(2));
Y = Y/(2*sqrt(2));
plot(X(:,1),X(:,2),'+');
hold on;
plot(Y(:,1),Y(:,2),'*');
epsilon = 1e-4;
gamma = log(5/epsilon)/3.7;

rmax = 2;
a = (exp(1)/(pi*Ngrid*Ncharges)^(1/4))*sqrt(gamma/2);

k = LogKernel(1);
A = Op(X,Y,k,a,1,epsilon);

% Fast way
fast = tic;
q_fast = A*V;
fast = toc(fast);

% Slow wat 
slow = tic;
q_slow = 0*X(:,1);
for i=1:Ncharges
    Yi_X1 = Y(i,1) - X(:,1);
    Yi_X2 = Y(i,2) - X(:,2);
    q_slow = q_slow + log(sqrt((Yi_X1).^2 + (Yi_X2).^2))*V(i);
end
slow = toc(slow);

err = q_slow - q_fast;
disp(max(abs(err)));

figure;
imagesc(Xaxis,Yaxis,real(reshape(q_fast,length(Xaxis),length(Yaxis))));

% Tracer une courbe de l'évolution de la complexité pour chaque méthode