clear all
close all


Ns = 30;

    
Ncharges = Ns;
Ngrid = 10000*Ncharges;

fprintf('%g charges on a %g grid\n\n',Ncharges,Ngrid);
fprintf('***************************************\n\n');
[X,Y,V,Xaxis,Yaxis] = GridAndUniformCircleData(Ngrid,Ncharges);
epsilon = 1e-4;
gamma = log(5/epsilon)/3.7;

rmax = 2;
a = min(0.05,(exp(1)/(pi*Ngrid*Ncharges)^(1/4))*sqrt(gamma/2));
k = nextY0root(30);
R1 = J0Kernel(k);
R2 = Y0Kernel(k);
assemble_i = tic;
A = -1i/4*(Op(X,Y,R1,a,1,epsilon) + 1i*Op(X,Y,R2,a,1,epsilon));
assemble = toc(assemble_i);


V = exp(1i*2*pi*rand(Ncharges,1)); %Random phases
V = V/norm(V,1);
q = A*V;



% q_val = 0*X(:,1);
% H=@(r)(-1i/4*(besselj(0,k*r) + 1i*bessely(0,k*r)));
% for j=1:Ncharges
%     Yi_X1 = Y(j,1) - X(:,1);
%     Yi_X2 = Y(j,2) - X(:,2);
%     q_val = q_val + H(sqrt((Yi_X1).^2 + (Yi_X2).^2))*V(j);
% end
% err = max(abs(q - q_val));
% fprintf('Maximal error: %g\n\n',err);

figure;
imagesc(Xaxis,Yaxis,(reshape(q.*conj(q),length(Xaxis),length(Yaxis))));
set(gca,'FontSize',14);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'dataAspectRatio',[1 1 1]);
% Tracer une courbe de l'évolution de la complexité pour chaque méthode

shadowTargetInds = find(X(:,1).^2 + X(:,2).^2 < 0.2); 



