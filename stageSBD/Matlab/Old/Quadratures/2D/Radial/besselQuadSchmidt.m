function[alpha,bessZs,fApprox,resL2] = besselQuadSchmidt(func,a,b,tol)

%% Applying Gram-Shcmidt method to get the quadrature

% Number of bessel modes computed in one loop
loadEvery = 10;


%% Initialization :

% Tolerance
tolq = tol^2;
reachedTol = false;
% Number of iterations
k = 0;
% k = q*loadEvery + t
q = 0;
t = 0;

% residual function
resFunc = @(X)(func(X));
alphaPrev = []; % func = sum(alphaPrev(k)*bessel(lam(k)*x)) + resFunc

% Orthonormolized Bessel coordinates (func(x) = sum(beta(k)*bessel_norm_k(x))
beta = zeros(loadEvery,1);

% Gram-Schmidt Matrix
B = []; % Coordinates of old vectors in the new basis
A = []; % A(i,j) = (ei|ej) in the new scalar product, ei being the old vectors
normSQ = []; % We store here the norms of the orthogonal system computed
% (for conditionning reasons).

%% Preparing loop

% Zeros of bessel function
bessZs = besselZeros(loadEvery);

% Quadrature to compute the scalar product of func against bessel functions
[x,w] = Gauss_Legendre1D(loadEvery*30,a,b);
ff = w.*x.*func(x);
% Square L2 norm of residue
resL2 = sum(w.*x.*resFunc(x).^2);
% Frequencies of bessel modes
lam = bessZs(q*loadEvery+1:(q+1)*loadEvery);
% quadrature scalar products of func against bessel modes
FF = sum((ones(size(lam))*ff').*besselj(0,lam*x'),2);


%% Main loop

while ~reachedTol
    
    % Updating results
    if t==loadEvery
        % Computing old bessel coordinates
        alpha = B*(normSQ*beta);
        % Old bessel coordinates of the residual
        dalpha = alpha - [alphaPrev;zeros(loadEvery,1)];
        % Updating alphaPrev
        alphaPrev = alpha;
        % Updating the quadrature of the residual
        resQuad = @(x)(sum(dalpha*ones(size(x(:)')).*...
            besselj(0,bessZs*x(:)'),1).');
        resFunc = @(x)(resFunc(x)-resQuad(x));
        
        % Changing quadrature resolution (more points)
        [x,w] = Gauss_Legendre1D((q+1)*loadEvery*30,a,b);     
        ff = w.*x.*func(x);
        resL2 = sum(w.*x.*resFunc(x).^2);
        
        q = q + 1;
        t = 0;
        
        bessZs = besselZeros((q+1)*loadEvery);
        lam = bessZs(q*loadEvery+1:(q+1)*loadEvery);
        FF = [FF; sum((ones(size(lam))*ff').*besselj(0,lam*x'),2)];
        beta = [beta; zeros(loadEvery,1)];
    end
    
    
    t = t+1;
    k = k+1;
    lamk = bessZs(k,1);
    
    % Computing A 
    for j = 1:k
        lamj = bessZs(j);
        A(j,k) = scalBess(lamj,lamk,a,b);
        A(k,j) = A(j,k);
    end
    % Next basis vector
    Ek = [ -(B*...
        (normSQ*...
        (normSQ*...
        (B'*...
        A(1:end-1,k)))));...
        1];
    
    % Norm of this vector
    normSQ(k,k) = 1/sqrt(Ek'*A*Ek);
    % Update B
    B = [[B;zeros(1,k-1)] Ek];
    % Coordinate along Ek
    beta(k,1) = normSQ(k,k)*sum(FF(1:k).*Ek);
    resL2 = resL2 - beta(k,1)^2;
    reachedTol = resL2 < tolq/2;
end


bessZs = bessZs(1:k);
beta = beta(1:k);
tt = linspace(a/100,1,1000);
figure
plot(tt,func(tt));
alpha = B*(normSQ*beta);
fApprox = @(x)(sum(alpha*ones(size(x(:)')).*...
    besselj(0,bessZs*x(:)'),1).');
hold on;
plot(tt,fApprox(tt));
title(sprintf('Radial approximation, %d components',k));



end
