function [alpha,rho,w0,reMakeError] = computeBesselCoeffH1_bis(a,b,func,derivative,tolInf,askGraph,Pmax,checkCond)


if ~exist('askGraph','var')
    askGraph = false;
end
if ~exist('Pmax','var')
    Pmax = Inf;
end
if ~exist('checkCond','var')
    checkCond = false;
end
if checkCond
    if tolInf ==0
        tolCond = 1;
    else
        tolCond = tolInf/10;
    end
end


fprintf('\n')

%% Defining parameters, initializing
timeTotal = tic;
gamma1_up = 0.3; gamma1_do = 0.28;
gamma2_up = -0.14; gamma2_do = -0.6;
gamma_up = max(gamma2_up+gamma1_up*log(1/tolInf),0.5);

UpperBound = min(fix(gamma_up*1/a)+1,Pmax);
LowerBound = min(max(fix((gamma2_do+gamma1_do*log(1/tolInf))*1/a)-1,1),Pmax);

fprintf('Number of components, first guess : %d \n',UpperBound);
firstFreqApprox = 3*pi/4; % Approximate of the next frequency to be computed


%% Enter in the main loop
% Pre-allocate rho :
% Compute next values of rho
tFreq = tic;
f  = @(r) besselj(0,r);
df = @(r) - besselj(1,r);
init_guess = firstFreqApprox + (0:UpperBound-1)'*pi;
rho = newton(f,df,init_guess);
tFreq = toc(tFreq);
fprintf(sprintf('%d frequencies computed in %5f seconds \n',UpperBound,tFreq));

%     % Check if only zeros
%
%     if max(abs(f(rho)))>1e-12
%         error('Problem in newton search for zeros')
%     end
%     % Check if all zeros
%     figure
%     t = linspace(0,rho(end),15*length(rho));
%     plot(t,besselj(0,t));
%     hold on
%     plot(rho,besselj(0,rho),'*r')
%     % Check complete
%


%% Store values to test for Linf norm.
tTestLinf = linspace(a,2*a,fix(100*UpperBound*a)+1);  % 50 points per freq
C = sqrt(2)./(abs(besselj(1,rho(:))));
D = besselj(0,rho(:)*[a b]);
E = besselj(1,rho(:)*[a b]);
storedValsJ0Mat = diag((C./rho))*besselj(0,rho(:)*tTestLinf(:)');
storedValsFunc = func(tTestLinf(:));

%% Compute A :
tic;
A = (-a*(C.*rho.*D(:,1))*(C.*E(:,1))'...
    + a*(C.*E(:,1))*(C.*rho.*D(:,1))'...
    +b*(C.*rho.*D(:,2))*(C.*E(:,2))'...
    -b*(C.*E(:,2))*(C.*rho.*D(:,2))')./(repmat(rho',UpperBound,1).^2-repmat(rho,1,UpperBound).^2);
A(1:(UpperBound+1):end) = C.^2.*(...
    b^2/2*(D(:,2).^2+E(:,2).^2) - b./rho.*D(:,2).*E(:,2) -...
    (a^2/2*(D(:,1).^2 + E(:,1).^2) - a./rho.*D(:,1).*E(:,1)));
computationA2 = toc;
fprintf('Matrix of scalar products computed in %f seconds \n',computationA2);

%% Compute the cholesky factorization
tChol = tic;
try
    fprintf('Cholesky : .... ')
    T = chol(A);
    fprintf('Done ! \n')
catch
    alpha = NaN;
    rho = NaN;
    w0 = NaN;
    return
end
% Check conditionning
if checkCond
    condCheck = cond(T)^2;
    if condCheck*eps > tolCond
        alpha = NaN;
        rho = NaN;
        w0 = NaN;
        return
    end
end

% B = inv(T')';
% disp(B(1,1));
tChol = toc(tChol);
fprintf('Inverse cholesky decomposition of A computed in %f seconds \n',tChol);

%% Compute the scalar products of func with basis functions
% Quadrature to compute the scalar product of func against bessel functions
tic
resolNGL = 5;
NGL = min(resolNGL*length(rho),4000);
[x,w] = Gauss_Legendre1D(NGL,a,b);
ff = w.*x.*derivative(x);

JJ = diag(C)*-besselj(1,rho(:)*x(:)');
f_ei = JJ*ff;
tBess = toc;

fprintf('Spent %f seconds evaluating Bessel functions \n',tBess)


beta = T'\f_ei;
storedValsOrthoMat = T'\storedValsJ0Mat;
% alphaTest = R\(Q.'*f_ei);
Beta = repmat(beta,1,UpperBound);
Beta = triu(Beta);
if nargout ==4
    reMakeError.beta = beta;
    reMakeError.T = T;
    reMakeError.Beta = Beta;
    reMakeError.storedValsOrthoMat = storedValsOrthoMat;
    reMakeError.storedValsFunc = storedValsFunc;
    reMakeError.A = A;
    
end

Pguess = LowerBound;
Beta = Beta(:,Pguess:end);
approx = storedValsOrthoMat'*Beta;
error = abs(repmat(storedValsFunc,1,size(approx,2))-approx);
error = max(error,[],1);
if isempty(find(error < tolInf, 1))
    P = Pmax;
else
    P = Pguess + find(error < tolInf,1,'first')-1;
end
alpha = (T(1:P,1:P)\beta(1:P))./rho(1:P);
w0 = alpha.*C(1:P);
timeTotal = toc(timeTotal);
rho = rho(1:P);

fprintf('Radial quadrature of %d components computed in %d seconds \n',P,timeTotal);

% Graph
fprintf('\n')
if askGraph
    figure
    quad = @(x)(besselApprox(rho,alpha,0,x,true));
    t = linspace(a/2,(b+1)/2,length(rho)*50);
    tTestLinf = linspace(a,b,length(rho)*50);
    plot(t,log(abs(func(t)-quad(t)+quad(b)-func(b))));
    figure
    plot(rho,alpha);
    
    resLinf = max(abs(quad(tTestLinf) - func(tTestLinf)));
    %     fprintf('Linf error : %d \n',resLinf);
    %
end



