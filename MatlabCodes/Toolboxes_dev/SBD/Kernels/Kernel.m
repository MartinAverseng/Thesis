classdef Kernel
    % Kernel
    
    properties (SetAccess = protected, GetAccess = public)
        func,der,scalFunc,normFunc,startFreq,gamma_est,customRadialQuad = [];
    end
    
    methods
        function[kernel] = Kernel(func,der)
            if nargin == 0
                func = @(x)(0*x);
                der = @(x)(0*x);
            end                
            kernel.func = func;
            kernel.der = der;
            fun1 = @(rho)(@(x)((x.*der(x)).*...
                -rho(:)'.*Cp(rho(:)').*besselj(1,rho(:)'*x)));
            % H10 scalar product
            fun2 = @(x)(x.*der(x).^2);
            % H10 norm
            kernel.scalFunc = @(a,b,rho)(2*pi*integral(fun1(rho),...
                a,b,'ArrayValued',true)); % -2\pi \int_{a}^b rf'(r)ep'(r)dr
            kernel.normFunc = @(a,b)(sqrt(2*pi*integral(fun2,a,b)));
            kernel.startFreq = 0; % We don't know a priori where the energy
            % is located in the spectrum. 
            kernel.gamma_est = @(tol)deal(0,7); % No fine tuning of gamma
            kernel.customRadialQuad = [];
        end
    end
    
    methods (Access = public)
        function[] = disp(this)
            fprintf('Kernel : function %s \n',func2str(this.func));            
        end
        function[this] = setScalFunc(this,f)
            this.scalFunc = f;
        end
        function[this] = setNormFunc(this,g)
            this.normFunc = g;
        end
        
        function[this] = setPBounds(this,h)
            this.PBounds = h;
        end
        function[this] = setStartFreq(this,k)
            this.startFreq = k;
        end
        function[C] = plus(k1,k2)
            C = Kernel(@(x)(k1.func(x) + k2.func(x)),@(x)(k1.der(x) + k2.der(x)));
            sf1 = k1.scalFunc;
            sf2 = k2.scalFunc;
            C.scalFunc =  @(a,b,rho)(sf1(a,b,rho) + sf2(a,b,rho));
            C.customRadialQuad = @(a,tol,varargin)(k1.radialQuadKernel(a,tol,varargin{:}) ...
                + k2.radialQuadKernel(a,tol,varargin{:}));
        end
        function[c] = mtimes(lambda,this)
            if and(isa(lambda,'double'),isscalar(lambda))
                assert(isa(this,'Kernel'));
                c = this;
                c.func = @(x)(lambda*this.func(x));
                c.der = @(x)(lambda*this.der(x));
                c.scalFunc = @(a,b,rho)(lambda*this.scalFunc(a,b,rho));
                c.normFunc = @(a,b)(abs(lambda)*this.normFunc(a,b));
                c.customRadialQuad = @(a,tol,varargin)(lambda*this.radialQuadKernel(a,tol,varargin{:}));
            else
                assert(isa(lambda,'double'),isscalar(lambda))
                c = times(this,lambda);
            end
        end
    end
    methods (Access = public)
        function[rq] = radialQuadKernel(this,a,tol,varargin)
            if ~ isempty(this.customRadialQuad)
                rq = this.customRadialQuad(a,tol,varargin{:});
            else
                rq = RadialQuadrature(a,this,tol,varargin{:});
            end
        end
        function[this] = dilatation(old,lambda)
            oldFunc = old.func;
            oldDer = old.der;
            ffunc = @(x)(oldFunc(lambda*x));
            dder = @(x)(lambda*oldDer(lambda*x));
            this = Kernel(ffunc,dder);
            this = this.setStartFreq(lambda*old.startFreq);
            this.customRadialQuad = @(a,tol,varargin)(dilatation(old.radialQuadKernel(a,tol,varargin{:}),lambda));
        end
        
        
    end
end

