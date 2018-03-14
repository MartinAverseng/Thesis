%% Filter frequency region

clear all;
close all;

Ncharges = 15;
Ngrid = 1000*Ncharges;

fprintf('%g charges on a %g grid\n\n',Ncharges,Ngrid);
fprintf('***************************************\n\n');
[X,Y,V,Xaxis,Yaxis] = GridAndUniformCircleData(Ngrid,Ncharges);

r = 0.1;
c = [0 0];
dist = (X(:,1)-c(1)).^2 + (X(:,2)-c(2)).^2;
shadowZone = find(dist< r^2);
dist = dist(shadowZone);
epsilon = 1e-4;
gamma = log(5/epsilon)/3.7;

a = min(0.03,(exp(1)/(pi*Ngrid*Ncharges)^(1/4))*sqrt(gamma/2));
f = 100;
omega = 2*pi*f;
c0 = 340;
k = omega/c0;
R1 = J0Kernel(k);
R2 = Y0Kernel(k);
assemble_i = tic;
Aplot = -1i/4*(Op(X,Y,R1,a,1,epsilon) + 1i*Op(X,Y,R2,a,1,epsilon));
Xshadow = X(shadowZone,:);
Aopt = -1i/4*(Op(Xshadow,Y,R1,a,1,epsilon) + 1i*Op(Xshadow,Y,R2,a,1,epsilon));
Aopttransp = 1i/4*(Op(Xshadow,Y,R1,a,1,epsilon).' - 1i*Op(Xshadow,Y,R2,a,1,epsilon).');
assemble = toc(assemble_i);

costFun = @(U)(costL2(U,Aopt,Aopttransp));
U0 = 2*pi*rand(Ncharges,1);
lb = zeros(Ncharges,1);
ub = ones(Ncharges,1);
options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
U = fmincon(costFun,U0,[],[],[],[],lb,ub,[],options);
costFun2 = @(U)(costLinf(U,Aopt,dist));
U = fmincon(costFun2,U,[],[],[],[],lb,ub);

V = exp(1i*2*pi*U);
V = V/norm(V,1);

q = Aplot*V;
power = 10*log(real(1/2*q.*conj(q)))/log(10);
power = power - max(power);
figure;
imagesc(Xaxis,Yaxis,reshape(power,length(Xaxis),length(Yaxis)));
set(gca,'dataAspectRatio',[1 1 1]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
t = linspace(0,2*pi,100);
X = c(1) + r*cos(t);
Y = c(2) + r*sin(t);
hold on;
plot(X,Y,'r');

V = exp(1i*2*pi*U0);
V = V/norm(V,1);
q = Aplot*V;
power = 10*log(real(1/2*q.*conj(q)))/log(10);
power = power - max(power);
figure;
imagesc(Xaxis,Yaxis,reshape(power,length(Xaxis),length(Yaxis)));
set(gca,'dataAspectRatio',[1 1 1]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
ylabel('Relative sound level (dB)')
