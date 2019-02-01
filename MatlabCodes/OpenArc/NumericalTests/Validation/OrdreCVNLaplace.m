%% Ordre de convergence simple couche segment Laplace
%% Create a mesh with non-uniform density and assemble the Somega

close all;
clear all;%#ok
clc;

x = @(t)(t);
y = @(t)(0*t);  
I = [-1,1];
curve = SimpleCurve(x,y,I);

k = 0;

Ns = [30 60 120 240 480 800];
n = 2;
lambda1 = R2toRfunc(@(Z)(4*Z(:,1).^2 - 1));
u01 = (n+1)/2*lambda1;
omegadxomega_u = -(n+1)*R2toRfunc(@(Z)(4*Z(:,1).^3 - 3*Z(:,1)));

[errH_12_nonreg,errU1,errU0] = testN_OrderCV(curve,lambda1, u01, Ns,k,'quadNum',10,'fullMatrix',true,'omdxomu',omegadxomega_u);
figure
% loglogTrislope(1./Ns(:),abs(errH_12_nonreg))
loglogTrislope(1./Ns(:),errU0(:))
hold on
loglogTrislope(1./Ns(:),1/10*errU1(:))

legend({'$L^2_\omega$ error','$ U^1 $ error'},'Interpreter','LaTeX');


set(gca,'XTick',1./(fliplr([50 100 200 400])))
xlabel('mesh size $h$','Interpreter','LaTeX')
ylabel('$L^2_\frac{1}{\omega}$ error $e_1(h)$','Interpreter','LaTeX')

set(gca,'XTick',1./(fliplr([50 100 200 400])))
xlabel('mesh size $h$','Interpreter','LaTeX')
ylabel('Approximation error','Interpreter','LaTeX')
legend('location','southeast');
set(gca,'Fontsize',15)