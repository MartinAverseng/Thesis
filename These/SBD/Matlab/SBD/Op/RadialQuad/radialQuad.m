function [out] = radialQuad(a,b,kernel,tol,varargin)
% This function computes the coefficients of the radial quadrature of a
% function which values and derivative are given by an object 'kernel' of
% type Kernel.
% A radial quadrature is an approximation of a radial function of the form 
% f(x) \approx \sum_{p=1}^P \alpha_p J_0(\rho_p |x|) 
% Valid for a < |x| < b at a tolerance tol for some coefficients \alpha_p and frequencies \rho_p. 
% This function applies a Gram-Schmidt method to orthonormalize the family
% of functions x -> J_0(\rho_p|x|) where \rho_p are the roots of J_0, 
% (this family is an orthonormal basis of H1_{0,rad}(B), B being the unit 
% ball) on the space H^1_{0,rad}(B(0,b) \B(0,a))
% No theory is available if b is ~= 1 (and often the method fails because
% of numerical instability in the Gram Matrix inversion) while for b=1,
% convergence at an exponential rate is ensured for the quadrature. 

if nargin ==0
    run('radialQuadVal.m')
    return
end
%% Timer initialization : 
tFreq = 0;
tStoreValForError = 0;
tAssembleA = 0;
tChol = 0;
tBess = 0;
timeCholInv = 0;
tTestErr = 0;
timeTotal = tic; % this one starts immediately :)

%% Arguments parsing
p = inputParser;
p.addRequired('a',@check_a);
p.addRequired('b',@(x)(check_b(x,a)));
p.addRequired('kernel',@(x)isa(x,'Kernel'));
p.addRequired('tol',@check_tol);

p.addOptional('Pmax',Inf,@(Pmax)(check_Pmax(Pmax)));
p.addOptional('verbose',false,@islogical);
p.addOptional('batch_size',1,@check_batch_size);

p.parse(a,b,kernel,tol,varargin{:});
vars = p.Results;
a = vars.a;
b = vars.b;
kernel = vars.kernel;
func = kernel.func;
derivative = kernel.der;
tol = vars.tol;
PBounds = vars.kernel.PBounds;
PBounds = PBounds(a,tol);
startFreq = vars.kernel.startFreq;
verbose = vars.verbose;
Pmax = vars.Pmax;

if b==1
    batch_size = fix(vars.batch_size*1/a)+1;
else
    batch_size = fix(vars.batch_size*(1/a))+1;
end

% Removing constant
if and(abs(func(b))>1e-12,~isa(kernel,'J0Kernel'))
    offset = func(b);
    func = @(x)(func(x) - func(b)); 
else
    offset = 0;
end

% Parse upper and lower bounds for number of components
PmaxLessThanHeuristicBound = false;
LowerBound = PBounds(1);
UpperBound = PBounds(2);
if UpperBound == Inf
    useHeuristicBounds = false;
    UpperBound = min(batch_size,Pmax);
else
    useHeuristicBounds = true;
    if Pmax < UpperBound
        PmaxLessThanHeuristicBound = true;
    end
    LowerBound = min(LowerBound,Pmax);
    UpperBound = min(UpperBound,Pmax);
end
if verbose&&(Pmax==Inf)
    fprintf('Number of components, first guess : %d \n',UpperBound);
end


% miscellaneous initialization
warnHeuristicNotWorked = false;
nIter = 0;
% Done parsing



%% Main Loop

if isa(kernel,'J0Kernel')
    % Very particuliar case of approximating J0 : no need to perform
    % projections !
    
    alpha0 = Inf; % I should never use this, if useful then make the effort of computing it 
    alpha = 1; % Here since func = J_0(Rx), thus trivial approximation !
    rho = kernel.R;
    resH10 = 0;
    errLinf = 0;
    reachedTol = true;
    scal01 = NaN; 
    nIter = 0;
else



satisfied = false;

while ~satisfied
    nIter = nIter + 1;
    assert(UpperBound <= Pmax);
    
    %% Compute frequencies
    
    tFreqTic = tic;
    rho = besselJroots(startFreq,UpperBound);
    tFreq = tFreq + toc(tFreqTic);
    
    if verbose
        fprintf(sprintf('%d frequencies computed in %5f seconds \n',UpperBound,tFreq));
    end
    
    %% Store values to test for Linf norm.
    
    tStoreValForErrorTic = tic;
    tTestLinf = linspace(a,min(2*a,b),fix(10*max(rho)*a)+1);  % 10 points per freq
    C = sqrt(2)./(abs(besselj(1,rho(:))));
    D = besselj(0,rho(:)*[a b]);
    E = besselj(1,rho(:)*[a b]);
    storedValsJ0Mat = diag((C./rho))*besselj(0,rho(:)*tTestLinf(:)');
    storedValsFunc = func(tTestLinf(:));
    tStoreValForError = tStoreValForError + toc(tStoreValForErrorTic);
    
    if verbose
        fprintf('Values of J0 stored to compute error in %f seconds \n',tStoreValForError);
    end
    
    %% Assemble matrix A :
    tAssembleATic = tic;
    A = (-a*(C.*rho.*D(:,1))*(C.*E(:,1))'...
        + a*(C.*E(:,1))*(C.*rho.*D(:,1))'...
        +b*(C.*rho.*D(:,2))*(C.*E(:,2))'...
        -b*(C.*E(:,2))*(C.*rho.*D(:,2))')./(repmat(rho',UpperBound,1).^2-repmat(rho,1,UpperBound).^2);
    A(1:(UpperBound+1):end) = C.^2.*(...
        b^2/2*(D(:,2).^2+E(:,2).^2) - b./rho.*D(:,2).*E(:,2) -...
        (a^2/2*(D(:,1).^2 + E(:,1).^2) - a./rho.*D(:,1).*E(:,1)));
    tAssembleA = tAssembleA + toc(tAssembleATic);
    if verbose
        fprintf('Matrix of scalar products computed in %f seconds \n',tAssembleA);
    end
    
    %% Compute Cholesky factorization of A
    tCholTic = tic;
    if verbose
        fprintf('Cholesky : .... ')
    end
    cholDone = false;
    while ~cholDone && batch_size>1
        try
            T = chol(A);
            cholDone = true;
        catch
            batch_sizeOld = batch_size;
            batch_size = round(batch_size/2);
            remove_size = batch_sizeOld - batch_size;
            A = A(1:end-remove_size,1:end-remove_size);
        end
    end
    if ~cholDone
        disp('Impossible to achieve required precision because of ill-conditionning of Gram Matrix')
        fprintf('Condition number : %s',num2str(cond(A)))
        error('Please decrease tolerance');
    end
    
    
    tChol = tChol + toc(tCholTic);
    if verbose
        fprintf('Done ! \n')
        fprintf('Inverse cholesky decomposition of A computed in %f seconds \n',tChol);
    end
    
    %% Compute H^1_0 scalar products on non-orthonormal basis
    
    tBessTic = tic;
    [f_ei,normH10,scal01] = scalProds(a,b,C,kernel,rho);
    tBess = tBess + toc(tBessTic);
    if verbose
        fprintf('Spent %f seconds evaluating scalar products \n',tBess)
    end
    
    %% Deduce orthonormal projections
    timeCholInvTic = tic;
    beta = T'\f_ei;
    % By definition of beta, 
    % f \approx \sum_{p=1}^UpperBound beta(p) \tilde{e}_p
    % where \tilde{e} is an orthonormal basis for the scalar product on
    % \mathcal{A}(a,b)
    
    
    
    %% Test for L^inf error and find P
    
    tTestErrTic = tic;
    storedValsOrthoMat = T'\storedValsJ0Mat;
    timeCholInv = timeCholInv + toc(timeCholInvTic);
    if verbose
        fprintf('Spent %f seconds solving for coefficients \n',timeCholInv)
    end
    Beta = repmat(beta,1,UpperBound);
    Beta = triu(Beta);
    Pguess = LowerBound;
    Beta = Beta(:,Pguess:end);
    approx = storedValsOrthoMat'*Beta;
    err = abs(repmat(storedValsFunc,1,size(approx,2))-approx);
    err = max(err,[],1);
    
    if isempty(find(err < tol, 1))
        % Tolerance not reached !
        reachedTol = false;
        
        if Pmax > UpperBound
            satisfied = false; % loop again 
            LowerBound = UpperBound;
            UpperBound = min(UpperBound + batch_size,Pmax);
        else
            satisfied = true; % we reached Pmax so we stop the computation
            P = Pmax; % = UpperBound
            errLinf = err(P - Pguess + 1);
        end
        
        if useHeuristicBounds && ~PmaxLessThanHeuristicBound
            warnHeuristicNotWorked = true;
        end
        
    else
        % Tolerance reached !

        reachedTol = true;
        satisfied = true;
        P = Pguess + find(err < tol,1,'first')-1;
        errLinf = err(P - Pguess + 1);
        assert(P <= Pmax);
    end
    tTestErr = tTestErr + toc(tTestErrTic);
    
    if verbose
        fprintf('Spent %f seconds evaluating error \n',tTestErr)
    end
    
end


%% Final coefficients (back to non-orthonormal basis)

beta = beta(1:P);rho = rho(1:P);T = T(1:P,1:P); C = C(1:P); scal01 = scal01(1:P);
assert(normH10^2 - sum(beta.^2) > -1e-10)
if normH10^2 > sum(beta.^2)
    resH10 = sqrt(normH10^2 - sum(beta.^2));
    if and(a>0,b==1)
        boundLinf = sqrt(-log(a)/(2*pi))*resH10;
        if errLinf > boundLinf
            error('Theroretical bound didn''t work. Either due to bad conditionning or error of code')
        end
    end
else
    resH10 = 0;
    warning(['Could not compute H10 residueal. Results may be inaccurate.'...
        ' If you entered a value of ''a'' above 0.3, consider decreasing ''a'' to improve stability'])
end
% Theoretical control : 


if warnHeuristicNotWorked
    warning('Heuristic bounds have been put in default')
    disp('Displaying arguments')
    disp(a); disp(b); disp(tol);disp(func);disp(derivative),disp(UpperBound);
    disp('End of Display')
end

if ~reachedTol
    warning('Tolerance not reached, because Pmax value was too small. Consider increasing Pmax')
end

%% Compute coefficients alpha

timeCholInvTic = tic;
alpha0 = (T\beta);
timeCholInv = timeCholInv + toc(timeCholInvTic);
alpha = alpha0.*(C./rho);

end

% Using alpha, the approximation takes the explicit form 
% f(x) = \sum_{p=1}^P alpha(p) besselj(0,rho(p)*x)
% Using alpha0, the approximation takes the form 
% f(x) = offset + \sum_{p=1}^P alpha0(p) e_p
% Where e_p are orthonormal in H10(B)

if verbose
    fprintf('Radial quadrature of %d components computed in %d seconds \n',P,timeTotal);
end

%% Timings
% Order important for consistence with timings display
timeTotal = toc(timeTotal);
times = struct;
times.freq = tFreq;
times.storeValsError = tStoreValForError;
times.assembleA = tAssembleA;
times.chol = tChol;
times.cholInv = timeCholInv;
times.projections = tBess;
times.testErr = tTestErr;
times.total = timeTotal;

%% Process output

out = struct;
out.a = a;
out.b = b;
out.func = func;
out.derivative = derivative;
out.tol = tol;
out.Pmax = Pmax;
out.startFreq = startFreq;
out.alpha0 = [offset; alpha0];
out.alpha = [offset; alpha];
out.rho = [0; rho];
out.resH10 = resH10;
out.errLinf = errLinf;
out.reachedTol = reachedTol;
out.times = times;
out.scal01 = scal01;
out.nIter = nIter;

end


%% Argument parsing functions
function[TF] = check_a(a)
if ~isscalar(a)
    error('Input is not scalar');
elseif ~isnumeric(a)
    error('Input is not numeric');
elseif (a < 0)
    error('Input must be >= 0');
elseif (a > 1)
    error('Input must be <= 1');
else
    TF = true;
end
end
function[TF] = check_b(b,a)
if ~isscalar(b)
    error('Input is not scalar');
elseif ~isnumeric(b)
    error('Input is not numeric');
elseif (b <= a)
    error('Input must be > a');
elseif (b > 1)
    error('Input must be <= 1');
else
    TF = true;
end
end
function[TF] = check_tol(tol)
if ~isscalar(tol)
    error('Input is not scalar');
elseif ~isnumeric(tol)
    error('Input is not numeric');
elseif (tol <= 0)
    error('Input must be > 0');
else
    TF = true;
end
end
function[TF] = check_Pmax(P)
if ~isscalar(P)
    error('Input is not scalar');
elseif ~or(mod(P,1)==0,P==Inf)
    error('Input is not  (integer / Inf)');
elseif (P <= 0)
    error('Input must be > 0');
elseif P ~= Inf
    TF = true;
end
end
function[TF] = check_batch_size(b)
if ~isscalar(b)
    error('Input is not scalar');
elseif ~(mod(b,1)==0)
    error('Input is not integer');
elseif (b <= 0)
    error('Input must be > 0');
else
    TF = true;
end
end


