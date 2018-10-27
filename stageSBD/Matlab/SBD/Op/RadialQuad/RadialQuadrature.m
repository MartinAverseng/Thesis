classdef RadialQuadrature
    % Approximation of a function as a series of the functions (e_i)
    % where e_i is the i-th normalized (in H10 norm) eigenfunction of the
    % Laplace operator with Dirichlet boundary conditions.
    % The approximation takes place on the set a < |x| < b of R^2
    % and writes
    % func \approx \sum_{i \in I} alpha(i) e_i(x)
    % The coefficients $\alpha$ are chosen as the minimizers of the H10
    % norm of this approximation.
    properties (GetAccess = public, SetAccess = protected)
        a, b, kernel, tol, Pmax, alpha0, alpha, rho,
        resH10, errLinf, reachedTol, times, scal01, nIter;
        
    end
    
    methods
        %% constructor and displays
        
        function[rq] = RadialQuadrature(aa,bb,kernel,ttol,varargin)
            %% Input description
            % REQUIRED INPUTS :
            % - a, b : as described in object description
            % - func : the radial function to be approximated
            % - der : the derivative of func
            % OPTIONAL INPUTS (Name-Value pairs)
            % - Pmax / type : integer or Inf / Default : Inf
            % Descr : maximal number of elements in quadrature
            % - verbose / type : logical / default : false
            % Descr : Set to true to enable text displays during computation
            % - batch_size : sets the constant K such that fix(K/a)+1 is the
            % number of frequencies added at each iteration.
            
            %% Creation of quadrature
           
            out = radialQuad(aa,bb,kernel,ttol,varargin{:});
            
            %% Export
            rq.a = out.a;
            rq.b = out.b;
            rq.kernel = kernel;
            rq.tol = out.tol;
            rq.Pmax = out.Pmax;
            rq.alpha0 = out.alpha0;
            rq.alpha = out.alpha;
            rq.rho = out.rho;
            rq.resH10 = out.resH10;
            rq.errLinf = out.errLinf;
            rq.reachedTol = out.reachedTol;
            rq.times = out.times;
            rq.scal01 = out.scal01;
            rq.nIter = out.nIter;
            
            rrho = rq.rho;
            assert(rrho(1) == 0); 
            % This convention must always be satisfied !! The constant. 
            
        end
        function[] = show(this)
            showAux(this);
        end
        function[] = disp(this)
            dispAux(this);
        end
        function[] = showTimes(this)
            t = this.times;
            fprintf('Time history (s) :\n\n')
            rowNames = {'J0 roots';'J0 values';'Assemble Gram';'Cholesky';'Cholesky^{-1}';'Projection';'Linf';'TOTAL'};
            Times = struct2array(t)';
            T = table(Times,'RowNames',rowNames,'VariableNames',{'t'});
            disp(T);
           
            explode = Times(1:end-1)*0+1;
             figure('Name','Computation time','NumberTitle','off')
            pie(Times(1:end-1)/sum(Times(1:end-1)),explode,rowNames(1:end-1))
        end
        
        %% Getters
        function[aa,bb] = getInterval(this)
            aa = this.a;
            bb = this.b;
        end
        function[ttol] = getTol(this)
            ttol = this.tol;
        end        
        function[s] = getQuadErrString(this)
            if this.errLinf < 0.1*this.tol
                s = sprintf('Quadrature error : < %s \n',num2str(this.errLinf));
            elseif this.errLinf < this.tol
                s= sprintf('Quadrature error : < %s \n',num2str(this.tol));
            else
                s = sprintf('Quadrature error : %s \n ',num2str(this.errLinf));
                s = [s,sprintf('The tolerance defined by user was not reached. Increase Pmax \n')];
            end
        end
        function[res] = getResH10(this)
            res = this.resH10;
        end
        function[res] = getTheoreticalBoundH10(this)
            res =  sqrt(-log(this.a)/(2*pi))*this.resH10;
        end
        function[alph1,rrho] = getAlpha(this)
            % Get the coefficients such that
            % func \approx \sum_{p=1}^P alph(p) besselj(0,rho(p)*x)
            alph1 = this.alpha;
            rrho = this.rho;
        end
        function[alph2,pVec] = getAlphaOrtho(this)
            % Get the coefficients alph2 and frequencies rrho such that
            % func \approx \sum_{p\in pVec} alph(p) e_p
            % where (e_p) is the orthonormal (in H10(B)) hilbert basis of
            % eigenvectors of the Laplace operator
            alph2 = this.alpha0;
            pVec = round(this.rho/pi+1/4);
            assert(isequal(pVec,unique(pVec)))
        end
        function[s] = getFunc(this)
            
                s = func2str(this.kernel.func);
        end        
        function[this] = dilatation(this,lambda)
            assert(lambda > this.a); % Otherwise new value of a would be > 1 
            this.a = this.a/lambda; 
            this.b = min(this.b/lambda,1);
            this.rho = this.rho*lambda;
            newKernel = this.kernel.dilatation(lambda);
            this.kernel = newKernel;
        end
    end
    
end

%% Auxiliary functions

function[] = dispAux(in)

func = in.kernel.func;
rho = in.rho;
a = in.a;
b = in.b;
P = length(rho);

if isa(in,'RadialQuadratureY0')
    s = 'Y0(Rx)';
    fprintf('Radial Quadrature of function %s \n',s);
    fprintf('With R = %s \n',num2str(in.getR));
elseif isa(in,'RadialQuadratureLog')
    s = 'log(Rx)';
    fprintf('Radial Quadrature of function %s \n',s);
    fprintf('With R = %s \n',num2str(in.getR));
else
    fprintf('Radial Quadrature of function %s \n',func2str(func));
end
fprintf('Domain of approximation : [%s, %s]\n',num2str(a),num2str(b))
fprintf('Number of components : %d \n',P)
disp(in.getQuadErrString);
in.showTimes();

end
function [] = showAux(in)

a = in.a;
b = in.b;
func = in.kernel.func;
rho = in.rho;
alpha = in.alpha;
alpha0 = in.alpha0;
tol = in.tol;
scal01 = in.scal01;
startFreq = in.kernel.startFreq;

crop = 0.99; % Needed to avoid having a 10^-16 value in the error at the outer edge
% Since function and approx are always equal at b
t2 = linspace(0,1,min((fix(max(rho)/pi)+1)*20,10000));
t1 = t2(and(t2>a/3,t2<b*crop));
quad1 = zeros(size(t1));
for p = 1:length(rho)
    quad1 = quad1 + alpha(p)*besselj(0,rho(p)*t1);
end
quad2 = zeros(size(t2));
for p = 1:length(rho)
    quad2 = quad2 + alpha(p)*besselj(0,rho(p)*t2);
end

figure('Name','Radial quadrature','NumberTitle','off')
subplot(1,3,1)
plot(t1,func(t1));
hold on
plot(t2,quad2,'--');
hold on
plot([a a],ylim,'k--');
if b<1
    plot([a a],ylim,'k--');
end
axis tight
title('Function and its radial quadrature','Interpreter','LaTex');
xlabel('r');

subplot(1,3,2)
loglog(t1,abs(func(t1)-quad1));
hold on
loglog(t1,t1*0+tol,'k--');
hold on
loglog([a a],ylim,'k--');
grid on
axis tight
title('Quadrature error','Interpreter','LaTex');
xlabel('r','Interpreter','LaTex');

subplot(1,3,3)
scatter(rho(2:end),abs(alpha0(2:end)).^2);
hold on

scatter(rho(2:end),abs(scal01).^2);

set(gca,'yscale','log');    

grid on;
axis tight
if startFreq >0
    hold on
    plot([startFreq startFreq],ylim,'k--','LineWidth',2);
end
title('Spectral power','Interpreter','LaTex');
if isempty(scal01(isnan(scal01)))
    legend({'Optimal approximation','Truncature of Bessel-Fourier expansion'})
end
xlabel('$\rho$','Interpreter','LaTex');

end


