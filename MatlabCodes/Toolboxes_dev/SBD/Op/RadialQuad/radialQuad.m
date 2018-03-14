function [rq] = radialQuad(a,kernel,tol,varargin)
% [rq] = radialQuad(a,kernel,tol,varargin)
% Computes the coefficients of the radial quadrature of a function which
% values and derivative are specified through an object 'kernel' of % type Kernel.
% The radial quadrature computed here is an approximation of the function of the form
% f(r) \approx \sum_p \alpha_p J_O(\rho_p r), a(1) < r < a(2)
% The coefficients \alpha_p are
% chosen as the minimizers of the H10 error of this approximation. The
% number of terms in the quadrature depends on the target error.
% We have shown that the approximation is sparse when there holds
% Delta f(1) = 0 (and even better if there exsists n>1 s.t. Delta^s f = 0 for all s <= n)
% This can be enforced by removing a constant times r^2 to the kernel.
% INPUTS
% - a : bounds within the approximation must be valid. 'a' can be specified
% either as a scalar (in which case a(2) is replaced by 1) or a size two
% vector. However, the linear system to be inverted is fully understood
% only when a(2) = 1. In general, if the approximation is needed on an interval
% [A,B], we suggest to rescale the kernel so that the approximation is sought on [a,1].
% You can later rescale the RadialQuadrature object accordingly. Ideally, the
% option to use a as a size-2 vector should never be needed.
% - kernel : a Kernel object, containing the fields 'func' and 'der'. If
% func(x) = log(x), use the LogKernel class, and if func(x) = Y_0(kx) use
% the Y0Kernel class, since they implement more efficient versions of the
% methods of the Kernel class.
% - tol : the required tolerance of the approximation. radialQuad finds the
% smallest P for which this tolerance is achieved.
% Optional arguments (Name-Value):
% - 'Pmax' : use Pmax to specify a maximal number of terms in the
% approximation. Default is Inf.
% - 'startFreq' : If the kernel's spectrum has energy near a specific
% frequency, we look for coefficients in this region in priority. Default is 0.
% - 'monitorDerivative' : If this is set to true, the algorithm stops as
% soon as the derivative is approximated up to the requested tolerance
% (instead of the values of the function itself).
% - 'verbose' : 0 (no text printed), 1 (some text printed), 2 (detailed
% output)
% OUTPUTS
% The output 'rq' contains the following fields :
% - a and b, the bounds of the interval of validity of the approximation
% - func and derivative (those of the kernel)
% - tol (the input)
% - Pmax (the input / default value)
% - rho : the frequencies of the approximation, s.t.
% f(r) \approx \sum_p \alpha_p J_O(rho(p) r)
% - startFreq (the input / default value)
% - alpha0 : the values of the coefficients such that
% f(r) \approx \sum_{p=1}^P alpha0(p) e_p(r)
% where e_p(r) = Cp(rho(p)) J_0(rho(p)r) (see functions Cp, and coeff2func).
% - alpha : the values of the coefficients such that
% f(r) \approx \sum_{p=1}^P alpha(p) J_0(rho(p) r)
% - resH10 : H10 norm of the residual error, that is
% sqrt{ 2\pi \int_{a(1)}^{a(2)} r(f(r) - \sum_{p=1}^P alpha0(p)e_p(r))^2}
% - errLinf : the maximal error of the approximation.
% - reachedTol : boolean indicating whether tolerance has been reached.
% - scal01 : values of 2\pi\int_{0}^{1} r(f'(r)ep'(r))dr
% - nIter : number of iterations the function has needed.

MAXIMAL_NUMBER_OF_COMPONENTS = 1500;
% You can change this value if you have enough memory.
% You can set it to Inf but it could lead to 'Out of Memory' errors.
if nargin ==0
    run('radialQuadVal.m')
    return
end


%% Timer initialization :
tFreq = 0;
tStoreValForError = 0;
tChol = 0;
tBess = 0;
timeCholInv = 0;
tTestErr = 0;
timeTotal = tic; % this one starts immediately :)

%% Arguments parsing
p = inputParser;
p.addRequired('a',@check_a);
p.addRequired('kernel',@(x)isa(x,'Kernel'));
p.addRequired('tol',@check_tol);

p.addOptional('monitorDerivative',false)
p.addOptional('Pmax',MAXIMAL_NUMBER_OF_COMPONENTS,@check_Pmax);
p.addOptional('startFreq',0,@check_startFreq);
p.addOptional('verbose',1);
p.KeepUnmatched = true;

p.parse(a,kernel,tol,varargin{:});
vars = p.Results;
a = vars.a;
a_stab = a;

if length(a)==2
    b = a(2);
    a = a(1);
else
    b = 1;
end

kernel = vars.kernel;
func = kernel.func;
der = kernel.der;
derivative = kernel.der;
tol = vars.tol;
startFreq = vars.startFreq;
monitorDerivative = vars.monitorDerivative;
monitorFunc = ~monitorDerivative;
verbose = vars.verbose;
Pmax = vars.Pmax;

if verbose>=1
    printStarLine;
    disp('SBD package : launching the radial decomposition')
end

if and(startFreq == 0,b==1)
    [lo,up] = kernel.gamma_est(tol);
    
    % Above this value, the gram matrix is ill-conditionned.
else
    % we have to guess.
    lo = 0;
    up = 7;
end

Pup = min(fix(up/a)+1,Pmax);
Plow = max(min(fix(lo/a),Pup-10),1);


% Removing constant
if and(abs(func(b))>1e-12,~isa(kernel,'J0Kernel'))
    offset = func(b);
    func = @(x)(func(x) - func(b));
else
    offset = 0;
end

nIter = 0;
% Done parsing


%% Main Loop
if a==1
    rho = [];
    % Cp is the normalization constant so that Cp*J0(rho |x|) has unit H10
    % norm (although not in H10 if R is not a root of J0, but anyway...)
    
    alpha0 = [];
    alpha = 1; % Here since func = J_0(Rx), trivial approximation !
    
    resH10 = 0;
    errLinf = 0;
    reachedTol = true;
    scal01 = NaN;
    nIter = 0;
elseif isa(kernel,'J0Kernel')
    if verbose>=1
        fprintf('This is a J0 Kernel : Bessel decomposition has 1 term \n')
    end
    % Very particuliar case of approximating J0 : no need to perform
    % projections !
    rho = kernel.R;
    C = Cp(rho);
    % Cp is the normalization constant so that Cp*J0(rho |x|) has unit H10
    % norm (although not in H10 if R is not a root of J0, but anyway...)
    
    alpha0 = 1./C;
    alpha = 1; % Here since func = J_0(Rx), trivial approximation !
    
    resH10 = 0;
    errLinf = 0;
    reachedTol = true;
    scal01 = NaN;
    nIter = 0;
else
    
    
    
    satisfied = false;
    
    while ~satisfied
        nIter = nIter + 1;
        assert(Pup <= Pmax);
        if verbose>=2
            fprintf('New iterations, with number of components : %s \n',num2str(Pup));
        end
        %% Compute frequencies
        
        tFreqTic = tic;
        if startFreq==0
            rho = besselJroots(0,Pup);
        else
            rho = unique([besselJroots(0,fix(Pup/2)+1);besselJroots(startFreq,fix(Pup/2)+1)],'stable');
            count = Pup;
            while length(rho) < Pup
                count = count*2;
                rho = unique([besselJroots(0,fix(count/2)+1);besselJroots(startFreq,fix(count/2)+1)],'stable');
            end
            Pup = length(rho);
        end
        
        tFreq = tFreq + toc(tFreqTic);
        
        if verbose>=2
            fprintf(sprintf('%d frequencies computed in %5f seconds \n',UpperBound,tFreq));
        end
        
        %% Store values to test for Linf norm.
        
        tStoreValForErrorTic = tic;
        tTestLinf = RadialQuadrature.tTestLinf(a,b,rho);
        if ~monitorFunc
            tTestLinf = tTestLinf(tTestLinf>a*1.05);
        end
        
        C = Cp(rho);
        if monitorFunc
            ep_test = repmat(C,1,length(tTestLinf)).*besselj(0,rho(:)*tTestLinf(:)');
            f_test = func(tTestLinf(:));
        else
            ep_test = -repmat(C.*rho,1,length(tTestLinf)).*besselj(1,rho(:)*tTestLinf(:)');
            f_test = der(tTestLinf(:));
        end
        
        tStoreValForError = tStoreValForError + toc(tStoreValForErrorTic);
        
        if verbose>=2
            fprintf('Values of J0 stored to compute error in %f seconds \n',tStoreValForError);
        end
        
        %% Get the gram Matrix cholesky factor.
        if and(a==0,b==1)
            T = eye(Pup);
        else
            [~,T,a_stab,tChol] = gramMatrix([a + (a==0)*eps,b],rho);
            if a_stab < a
                Plow = 1; % We retry all values.
            end
        end
        %% Compute H^1_0 scalar products on non-orthonormal basis
        
        tBessTic = tic;
        [f_ei,normH10,scal01] = scalProds([a_stab,b],kernel,rho);
        tBess = tBess + toc(tBessTic);
        if verbose>=2
            fprintf('Spent %f seconds evaluating scalar products \n',tBess)
        end
        
        %% Deduce orthonormal projections
        timeCholInvTic = tic;
        beta = T'\f_ei;
        % By definition of beta,
        % f \approx \sum_{p=1}^P beta(p) \tilde{e}_p
        % where \tilde{e} is an orthonormal basis for the scalar product on
        % \mathcal{A}(a,b)
        
        
        
        %% Test for L^inf error and find P
        
        tTestErrTic = tic;
        storedValsOrthoMat = T'\ep_test;
        timeCholInv = timeCholInv + toc(timeCholInvTic);
        if verbose>=2
            fprintf('Spent %f seconds solving for coefficients \n',timeCholInv)
        end
        Beta = repmat(beta,1,Pup);
        Beta = triu(Beta);
        Beta = Beta(:,Plow:end);
        approx = storedValsOrthoMat'*Beta;
        err = abs(repmat(f_test,1,size(approx,2))-approx);
        err = max(err,[],1);
        
        if isempty(find(err < tol, 1))
            % Tolerance not reached !
            reachedTol = false;
            
            if Pmax > Pup
                satisfied = false; % loop again
                Plow = Pup; % -10 Just in case 'a' changes.
                Pup = min(2*Pup,Pmax);
            else
                satisfied = true; % we reached Pmax so we stop the computation
                Pup = Pmax;
                errLinf = err(Pup - Plow + 1);
            end
            
            
        else
            % Tolerance reached !
            
            reachedTol = true;
            satisfied = true;
            Pup = Plow + find(err < tol,1,'first')-1;
            errLinf = err(Pup - Plow + 1);
            assert(Pup <= Pmax);
            if verbose >=1
                fprintf('Bessel decomposition successfully computed. Number of terms : %s \n',num2str(Pup))
            end
        end            
        tTestErr = tTestErr + toc(tTestErrTic);
        
        if verbose>=2
            fprintf('Spent %f seconds evaluating error \n',tTestErr)
        end
        
    end
    
    
    %% Final coefficients (back to non-orthonormal basis)
    
    beta = beta(1:Pup);rho = rho(1:Pup);T = T(1:Pup,1:Pup); C = C(1:Pup); scal01 = scal01(1:Pup);
    %assert(normH10^2 - sum(beta.^2) > -1e-10)
    if normH10^2 > sum(beta.^2)
        resH10 = sqrt(normH10^2 - sum(beta.^2));
        if and(a>0,b==1)
            if monitorFunc
                boundLinf = sqrt(-log(a)/(2*pi))*resH10;
                if errLinf > boundLinf
                    warning('Theroretical bound didn''t work. Either due to bad conditionning or error of code')
                    assert(~ any(isa(kernel,'LogKernel'),isa(kernel,'Y0Kernel')),'This should never happen in the case of standard kernels.')
                end
            else
                
            end
        end
    else
        resH10 = 0;
        warning(['Could not compute H10 residueal. Results may be inaccurate.'...
            'This might be due to parameter ''a'' set too high.'])
    end
    if ~reachedTol
        assert(Pup==Pmax)
        if Pmax == MAXIMAL_NUMBER_OF_COMPONENTS
            warning(['The DEFAULT max number of components is too small to '...
                'reach requested tolerance']);
        else
            warning('Tolerance not reached, because Pmax value was too small. Consider increasing Pmax, or decreasing a')
        end
    end
    
    %% Compute coefficients alpha
    
    timeCholInvTic = tic;
    alpha0 = (T\beta);
    timeCholInv = timeCholInv + toc(timeCholInvTic);
    alpha = alpha0.*C;
    
end

% Using alpha, the approximation takes the explicit form
% f(x) = \sum_{p=1}^P alpha(p) besselj(0,rho(p)*x)
% Using alpha0, the approximation takes the form
% f(x) = offset + \sum_{p=1}^P alpha0(p) e_p
% Where e_p are orthonormal in H10(B)



%% Timings
% Order important for consistence with timings display
timeTotal = toc(timeTotal);
if verbose >=1
    fprintf('Radial quadrature of %s components computed in %s seconds \n',num2str(Pup),num2str(timeTotal));
end
times = struct;
times.freq = tFreq;
times.storeValsError = tStoreValForError;
times.chol = tChol;
times.cholInv = timeCholInv;
times.projections = tBess;
times.testErr = tTestErr;
times.total = timeTotal;

%% Process output

rq = struct;
rq.a = a;
rq.b = b;
rq.func = func;
rq.derivative = derivative;
rq.tol = tol;
rq.Pmax = Pmax;
rq.startFreq = startFreq;
rq.alpha0 = [offset; alpha0];
rq.alpha = [offset; alpha];
rq.rho = [0; rho];
rq.resH10 = resH10;
rq.errLinf = errLinf;
rq.reachedTol = reachedTol;
rq.times = times;
rq.scal01 = scal01;
rq.nIter = nIter;
rq.verbose = verbose;
if verbose >=1
    fprintf('Done \n\n')
end
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
function[TF] = check_tol(tol)
if ~isscalar(tol)
    error('Input is not scalar');
elseif ~isnumeric(tol)
    error('Input is not numeric');
else
    TF = true;
end
end
function[TF] = check_startFreq(sf)
if ~isscalar(sf)
    error('Input is not scalar');
elseif ~isnumeric(sf)
    error('Input is not numeric');
elseif (sf <= 0)
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


