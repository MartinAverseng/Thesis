classdef LogKernel < Kernel
    % Object of type kernel but optimized for log so that the computation of the 
    % radial quadrature goes faster
    properties (Access = public)
        R;  %such that kernel =  x -> log(Rx) 
    end
    
    methods
      	function[kernel] = LogKernel(RR)
            kernel@Kernel(@(x)(log(RR*x)),@(x)(1./x))
            kernel = kernel.setScalFunc(@LogKernel.scalFuncLaplace);
            kernel = kernel.setNormFunc(@LogKernel.normFuncLaplace);
            kernel.gamma_est = @LogKernel.gamma_est;
            kernel.R = RR;
        end
    end
    methods (Access = public)
        function[this] = dilatation(this,lambda)
            this = LogKernel(lambda*this.R);    
            % No changes in gamma_est, scalFunc nor normFunc
        end 
    end
    methods (Access = protected)
         
    end
    methods (Static, Access = protected)
        function[low,up] = gamma_est(tol)
            % Helps the radial quadrature to guess the number of components
            up = 1/5.8*log(180/tol);
            low= max(up-1,0.5);
        end
        function[res] = scalFuncLaplace(a,b,rho)
            res = laplaceSP([a,b],rho);            
        end
        function[res] = normFuncLaplace(a,b)
            % Computes (a,b) -> sqrt(2*pi*\int_a^b x [(log(x))']^2 dx))
            % (= sqrt(2*pi*\int_a^b 1/x dx) = log(b/a)
            res = sqrt(2*pi*log(b/a));
        end
              
    end
    
end

