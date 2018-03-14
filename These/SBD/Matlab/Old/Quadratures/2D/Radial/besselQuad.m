function[alpha,bessZs,fApprox] = besselQuadSchmidt(func,a,b,tol)

func2 = @(X)(func(X));
loadEvery = 20;
bessZs = besselZeros(loadEvery);
[x,w] = Gauss_Legendre1D(loadEvery*2,a,b);
reachedTol = false;
k = 0;
q = 0;
t = 0;
ff = w.*x.*func2(x);
lam = bessZs(q*loadEvery+1:(q+1)*loadEvery);
FF = sum((ones(size(lam))*ff').*besselj(0,lam*x'),2);
beta = zeros(loadEvery,1);
B = [];
A = [];


while ~reachedTol
    k = k+1;
    t = t+1;
    if t==loadEvery
        q = q + 1;
        t = 0;
        bessZs = besselZeros((q+1)*loadEvery);
        lam = bessZs(q*loadEvery+1:(q+1)*loadEvery);
        [x,w] = Gauss_Legendre1D((q+1)*loadEvery*2,a,b);
        ff = w.*x.*func2(x);
        
        FF = [FF; sum((ones(size(lam))*ff').*besselj(0,lam*x'),2)];
        beta = [beta; zeros(loadEvery,1)];
    end
    lamk = bessZs(k,1);
    for j = 1:k
        lamj = bessZs(j);
        A(j,k) = scalBess(lamj,lamk,a,b);
        A(k,j) = A(j,k);
    end
    Ek = [-B*(B'*A(1:end-1,k));1];
    
    normSQ = Ek'*A*Ek;        
    Ek = Ek/sqrt(normSQ);
    B = [[B;zeros(1,k-1)] Ek];
    beta(k,1) = sum(FF(1:k).*Ek);
    reachedTol = abs(beta(k)) < tol/5;
    
    
end
bessZs = bessZs(1:k);
beta = beta(1:k);
tt = linspace(a/100,1,1000);
figure
plot(tt,func(tt));
alpha = B*beta;
fApprox = @(x)(sum(alpha*ones(size(x(:)')).*...
    besselj(0,bessZs*x(:)'),1).');
hold on;
plot(tt,fApprox(tt));
title(sprintf('Radial approximation, %d components',k));



end
