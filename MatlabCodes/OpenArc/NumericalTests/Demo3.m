%% Testing the preconditioners on non flat arcs and non-zero frequency
% We test the new preconditioning technique on a spirale

clear all %#ok
close all
clc;
k = 13;
arc = spirale;
figure 
subplot(1,2,1);
plot(arc)
title('The Spirale geometry');

N = 800;
repartition = @cos;
bounds = [-pi,0];
mesh = MeshCurve(arc,N,repartition,bounds);
subplot(1,2,2);
plot(mesh);
title('The adapted mesh');

Vh = weightedFEspace(mesh,'P1','1/sqrt(1-t^2)',5);
M = full(Vh.Mass);
Wh =  weightedFEspace(mesh,'P1','sqrt(1-t^2)',5);
dM = full(Wh.dMass);
omega2 = full(Wh.Mass);
%% 
% We assemble the single layer on this space and the three preconditioners :

Somega = singleLayer(k,Vh,[],{'full',true});
Sgalerk = Somega.galerkine(Vh,'U');
Mat = 4*M^(-1)*full(Sgalerk)*M^(-1)*(dM-k^2*omega2)*M^(-1);
Mat_pourrie = 4*M^(-1)*full(Sgalerk)*M^(-1)*(dM+1)*M^(-1);
EV = eig(Mat*full(Sgalerk));
C = cond(Mat*full(Sgalerk))
C = cond(Mat_pourrie*full(Sgalerk))
figure
plot(real(EV),imag(EV),'*');

Prec1 = Sgalerk.concretePart;
Prec2 = @(u)(M\(Sgalerk*(M\(dM*(M\u)))));
Prec3 = @(u)(M\(Sgalerk*(M\((dM-k^2*omega2)*(M\u)))));
Z0 = [-0.9,0.01];
u0 = R2toRfunc(@(Z)(log(sqrt((Z(:,1)-Z0(:,1)).^2 + (Z(:,2)-Z0(:,2)).^2)).*sin(3*k*Z(:,1))));
l = Vh.secondMember(u0);
fprintf('\nWithout Preconditioner\n');
tic;
[lambda0,~,~,~,resvec0]  = variationalSol(Sgalerk,l,40,1e-8,N);
fprintf('\n%s seconds ! \n\n S0\n',num2str(toc));
tic;
[lambda1,~,~,~,resvec1]  = variationalSol(Sgalerk,l,40,1e-8,N,Prec1);
fprintf('\n%s seconds !\n\n -SomegaDelta\n',num2str(toc));
tic;
[lambda2,~,~,~,resvec2]  = variationalSol(Sgalerk,l,[],1e-8,N,Prec2);
fprintf('\n%s seconds !\n\n -Somega(Delta+k^2) \n',num2str(toc));
tic;
[lambda3,~,~,~,resvec3]  = variationalSol(Sgalerk,l,[],1e-8,N,Prec3);
fprintf('\n%s seconds !\n',num2str(toc));
figure
semilogy(0:length(resvec0)-1,resvec0/norm(full(l)),'DisplayName','No precond')
hold on
semilogy(0:length(resvec1)-1,resvec1/norm(Prec1\full(l)),'DisplayName','S_0')
semilogy(0:length(resvec2)-1,resvec2/norm(Prec2(full(l))),'DisplayName','-S_{\omega}\Delta_{\omega}')
semilogy(0:length(resvec3)-1,resvec3/norm(Prec3(full(l))),'DisplayName','-S_{\omega}(\Delta_{\omega}+ k^2)')
legend show


