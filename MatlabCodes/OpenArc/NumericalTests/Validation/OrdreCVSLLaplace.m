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

Ns = [30 60 120 240 480];

X = R2toRfunc(@(Z)(Z(:,1)));

lambda1 = sqrt(1 - X^2);
u01 = -1/(2*pi)*((1 - X)*log(1 - X) + (1 + X)*log(1 + X) - 2);

[errH_12_nonreg,errL2_nonreg] = testSL_OrderCV(curve,lambda1, u01, Ns,k,...
    'correcMethod','constantTerm','quadNum',15,'fullMatrix',true);
figure
% loglogTrislope(1./Ns(:),abs(errH_12_nonreg))
% hold on
loglogTrislope(1./Ns(:),errL2_nonreg)

lambda2 = (1 - X^2)^(3/2);
u02 = 2/3*log(1 + X) - 1/3*X^3*log(X + 1) - 16/9 + X*log(X + 1) + 2/3*X^2 + 2/3*log(1 - X)+1/3*X^3*log(1 - X) - X*log(1 - X);
u02 = -1/(2*pi)*(u02);

[errH_12,errL2] = testSL_OrderCV(curve,lambda2, u02, Ns,k,...
     'correcMethod','constantTerm','quadNum',15,'fullMatrix',true);
% figure
% % loglogTrislope(1./Ns(:),real(errH_12))
hold on
loglogTrislope(1./Ns(:),errL2)

legend({'$\alpha = \alpha_1 = \omega \in T^{3/2}$','$\alpha = \alpha_2 = \omega^3 \in T^2$'},'Interpreter','LaTeX');
legend('location','southeast')
set(gca,'Fontsize',15)

set(gca,'XTick',1./(fliplr([50 100 200 400])))
xlabel('mesh size $h$','Interpreter','LaTeX')
ylabel('$L^2_\frac{1}{\omega}$ error $e_1(h)$','Interpreter','LaTeX')