clear all
close all


Ns = [15 50 150 500 1500 5000 15000];

for i = 1:length(Ns)
    
Ncharges = Ns(i);
Ngrid = 100*Ncharges;

fprintf('%g charges on a %g grid\n\n',Ncharges,Ngrid);
fprintf('***************************************\n\n');
[X,Y,V,Xaxis,Yaxis] = GridAndCloud(Ngrid,Ncharges);
X = X/(2*sqrt(2));
Y = Y/(2*sqrt(2));
epsilon = 1e-4;
gamma = log(5/epsilon)/3.7;

rmax = 2;
a = (exp(1)/(pi*Ngrid*Ncharges)^(1/4))*sqrt(gamma/2);

k = LogKernel(1);
assemble_i = tic;
A = Op(X,Y,k,a,1,epsilon);
assemble(i) = toc(assemble_i);

fprintf('Operator assembled in %g s\n\n',assemble(i))
% Fast way
fast_i = tic;
q_fast = A*V;
fast(i) = toc(fast_i);
fprintf('Fast product computed in %g s\n\n',fast(i))
% Slow wat 
slow_i = tic;
q_slow = 0*X(:,1);
for j=1:Ncharges
    Yi_X1 = Y(j,1) - X(:,1);
    Yi_X2 = Y(j,2) - X(:,2);
    q_slow = q_slow + log(sqrt((Yi_X1).^2 + (Yi_X2).^2))*V(j);
end
slow(i) = toc(slow_i);

fprintf('Brute force interactions computed in %g s\n\n',slow(i))
err(i) = max(abs(q_slow - q_fast));
fprintf('Maximal error : %g \n\n',err(i));


end

figure;
imagesc(Xaxis,Yaxis,real(reshape(q_fast,length(Xaxis),length(Yaxis))));

% Tracer une courbe de l'évolution de la complexité pour chaque méthode