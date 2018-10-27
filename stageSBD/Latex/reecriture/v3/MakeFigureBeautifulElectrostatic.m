clear all
close all


Ns = 50000;

    
Ncharges = Ns;
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
assemble = toc(assemble_i);


V = randn(Ns,1);
q = A*V;
figure;
imagesc(Xaxis,Yaxis,real(reshape(q,length(Xaxis),length(Yaxis))));
set(gca,'FontSize',14);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'dataAspectRatio',[1 1 1]);
% Tracer une courbe de l'évolution de la complexité pour chaque méthode