%% Test 1

clear all;
clc
close all;

Nx = 100;
X1 = (1:Nx)';
X2 = round(rand(Nx,1));
X = [X1,X2];
Y = [X1,1-X2];
V = randn(Nx,1);
V = V/norm(V,1);
dist = @(x,x0) max( sqrt( (x(:,1)-x0(1)).^2 + (x(:,2)-x0(2)).^2) );
rX   = dist(X,mean(X));
rY   = dist(Y,mean(Y));
rXY  = dist(mean(X),mean(Y));
rMax = rXY + rX  + rY;
rMin = 1;
x = X/rMax;
y = Y/rMax;
a = rMin/(rMax);
tol = 1e-3;
kernel = LogKernel(rMax);
rq = kernel.radialQuadKernel(a,tol);
q2D = Quad2D(rq);
q = q2D.conv(x,y,V);

% Validation 
qVal = zeros(Nx,1);
for i = 1:Nx
    for j = 1:Nx
        qVal(i) = qVal(i) + log(norm(X(i,:)-Y(j,:),2))*V(j);
    end
end
figure
plot(real(q));
hold on
plot(qVal,'--');
assert(norm(qVal - q,'inf')<=tol);


%% Test 2 

clear all;
clc
close all;

Nx = 150;
X1 = (1:Nx)';
X2 = round(rand(Nx,1));
X = [X1,X2];
Y = [X1,1-X2];
V = randn(Nx,1);
V = V/norm(V,1);
dist = @(x,x0) max( sqrt( (x(:,1)-x0(1)).^2 + (x(:,2)-x0(2)).^2) );
rX   = dist(X,mean(X));
rY   = dist(Y,mean(Y));
rXY  = dist(mean(X),mean(Y));
rMax = rXY + rX  + rY;
rMin = 1;
x = X/rMax;
y = Y/rMax;

a = rMin / rMax;
tol = 1e-3;

k = 10;
kernel = Y0Kernel(k*rMax);
rq = kernel.radialQuadKernel(a,tol);
q2D = Quad2D(rq);
q = q2D.conv(x,y,V);

% Validation 
for i = 1:Nx
    qVal(i,1) = 0;
    for j = 1:Nx
        qVal(i) = qVal(i) + bessely(0,norm(k*X(i,:)-k*Y(j,:),2))*V(j);
    end
end
figure
plot(real(q));
hold on
plot(qVal,'--');
assert(norm(qVal - q,'inf')<=tol);

%% 
Main;
a = 0.012;
tol = 1e-3;
kernel = Y0Kernel(400);
rq = kernel.radialQuadKernel(a,0,'Pmax',fix(5/a));
q2D = Quad2D(rq);
[dx,dy] = grad(q2D);

r = [linspace(a/2,a-0.001,500) linspace(a,1,1000)]';
figure;

for theta = 0:0.1:2*pi
    x = r*cos(theta); y = r*sin(theta);
    u = [x y];
    true = kernel.func(r);
    true_dx = kernel.der(r).*x./r;
    true_dy = kernel.der(r).*y./r;
    approx = q2D.value(u);
    approx_dx = dx.value(u);
    approx_dy = dy.value(u);
    subplot(1,3,1)
    loglog(r,abs(true_dx - approx_dx))
    subplot(1,3,2)
    loglog(r,abs(true_dy - approx_dy))
    subplot(1,3,3)
    loglog(r,abs(true - approx))
    pause(2)
end




%% 
clc 
close all
disp('success')