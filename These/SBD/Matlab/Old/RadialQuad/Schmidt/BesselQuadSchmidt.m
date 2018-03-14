function [ alpha,bessZs,fApprox,resL2,beta,B,normSQ,timeQuad ] = BesselQuadSchmidt( D,c,func,a,b,tol,freqCenter,Pmax )
% Computes the coefficients of the Fourier-Bessel Series
% giving the approximation 'quad' of 'func' satisfying the following criteria :
% - ||func - quad|| < eps
% - quad as the minimal number of coefficients
% The coefficients used in the series are the eigenvalues of the laplacians
% in dimension D, obeying to a Robin condition of the form  f'(x) + cf(x) = 0
% on the boundary of the unit ball. If c is set to infty, the condition
% above is replaced by Dirichlet condition.


if ~exist('freqCenter','var')
    freqCenter =0;
end
if ~exist('Pmax','var')
    Pmax =Inf;
end

nu = D/2 - 1;
tic;


%% Initialization :

loadEvery = 100;
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
Aprime = [];
lambdaMin = [];
normSQ = []; % We store here the norms of the orthogonal system computed
% (for conditionning reasons).

%% Preparing loop

% Zeros of bessel function
bessZs = BesselZeros(c,nu,loadEvery,freqCenter);



% Quadrature to compute the scalar product of func against bessel functions
[x,w] = Gauss_Legendre1D(4000,a,b);
ff = w.*x.^(D/2).*func(x);
% Square L2 norm of residue
resL2 = sum(w.*x.^(D-1).*resFunc(x).^2);
resI = resL2;
% Frequencies of bessel modes
lam = bessZs(q*loadEvery+1:(q+1)*loadEvery);
% quadrature scalar products of func against bessel modes
FF = sum((ones(size(lam))*ff').*(besselj(nu,lam*x')./(lam*ones(size(x'))).^nu),2);


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
            besselj(nu,bessZs*x(:)')./(bessZs*x(:)').^nu,1).');
        resFunc = @(x)(resFunc(x)-resQuad(x));
        
        % Changing quadrature resolution (more points)
%         [x,w] = Gauss_Legendre1D(2000,a,b);
        ff = w.*x.^(D/2).*func(x);
        resL2 = sum(w.*x.^(D-1).*resFunc(x).^2);
        
        q = q + 1;
        t = 0;
        
        bessZs = BesselZeros(c,nu,(q+1)*loadEvery,freqCenter);
        lam = bessZs(q*loadEvery+1:(q+1)*loadEvery);
        FF = [FF; sum((ones(size(lam))*ff').*(besselj(nu,lam*x')./(lam*ones(size(x'))).^nu),2)];
        beta = [beta; zeros(loadEvery,1)];
    end
    
    
    t = t+1;
    k = k+1;
    lamk = bessZs(k,1);
    
    % Computing A
    for j = 1:k
        lamj = bessZs(j);
        A(j,k) = ScalBess(lamj,lamk,nu,a,b)/(lamj*lamk)^nu;
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
    % Computing old bessel coordinates
        alpha = B*(normSQ*beta(1:k));
        % Old bessel coordinates of the residual
        % Updating alphaPrev
        % Updating the quadrature of the residual
        quad = @(x)(besselApprox(bessZs(1:k),alpha,0,x,0));
        resLinf = abs(quad(2*a)-func(2*a))
        reachedTol = resLinf < tol;
 end

timeQuad = toc;
bessZs = bessZs(1:k);
beta = beta(1:k);
tt = linspace(a/2,1,20);
figure
plot(tt,func(tt));
alpha = B*(normSQ*beta);
fApprox = @(x)(sum(alpha*ones(size(x(:)')).*...
    besselj(nu,bessZs*x(:)')./(bessZs*x(:)').^nu,1).');
hold on;

tt = linspace(0,1,1000);
plot(tt,fApprox(tt));
title(sprintf('Radial approximation, %d components',k));

figure
scatter(bessZs,alpha);
title('Coefficients of the decomposition')
xlabel('Frequency')
ylabel('Amplitude')

figure
scatter(1:length(bessZs),abs(alpha*sqrt(2)).^(-2/3)./(1:length(bessZs))'/pi);
title('Inverse of the oefficients in the decomposition')
xlabel('Frequency')
ylabel('Amplitude (dB)')

figure
plot(tt,log(abs(fApprox(tt)'-func(tt))));
title('Log error of the decomposition')
xlabel('x')
ylabel('error (dB)')

P = length(alpha); % number of components in quadrature
fprintf('Radial quadrature : %d components, \n %s seconds \n %d L2 error  \n\n',P,num2str(timeQuad),sqrt(resL2));

end

