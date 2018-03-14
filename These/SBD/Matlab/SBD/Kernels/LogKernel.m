classdef LogKernel < Kernel
    % Object of type kernel but optimized for log so that the computation of the 
    % radial quadrature goes faster
    properties (Access = public)
        R;  %such that kernel =  x -> Y0(Rx) 
    end
    
    methods
      	function[kernel] = LogKernel(RR)
            kernel@Kernel(@(x)(log(RR*x)),@(x)(1./x))
            kernel = kernel.setScalFunc(@LogKernel.scalFuncLaplace);
            kernel = kernel.setNormFunc(@LogKernel.normFuncLaplace);
            kernel = kernel.setPBounds(@LogKernel.PBounds);
            kernel.R = RR;
        end
    end
    methods (Access = public)
        function[this] = dilatation(this,lambda)
            this = LogKernel(lambda*this.R);            
        end 
    end
    methods (Access = protected)
         
    end
    methods (Static, Access = protected)
        function[loup] = PBounds(a,tol)
            % Helps the radial quadrature to guess the number of components
            loup = HeuristicBounds(a,tol);
        end
        function[res] = scalFuncLaplace(a,b,rho)
            % Computes (a,b,\rho) -> \int_a^b x(log(x))'J_1(\rho x) dx
            % ( = \int_a^b J_1(\rho x) dx )
            res = (besselj(0,rho*a) - besselj(0,rho*b))./rho;
            res = res(:)';
            % Explicit form :
            % https://www.wolframalpha.com/input/?i=integral+from+a+to+b+of++besselj(1,rho+x)
        end
        function[res] = normFuncLaplace(a,b)
            % Computes (a,b) -> \int_a^b x [(log(x))']^2 dx)
            % (= \int_a^b 1/x dx = log(b/a)
            res = sqrt(log(b/a));
        end
              
    end
    
end

