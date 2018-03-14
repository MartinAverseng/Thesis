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
b = 1;
tol = 1e-3;
kernel = LogKernel(rMax);
q2D = kernel.quad2D(a,b,tol);
q = q2D.conv(x,y,V);

% Validation 
for i = 1:Nx
    qVal(i,1) = 0;
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
b = 1;
tol = 1e-3;

k = 10;
kernel = Y0Kernel(k*rMax);
q2D = kernel.quad2D(a,b,tol);
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
clc 
close all
disp('success')