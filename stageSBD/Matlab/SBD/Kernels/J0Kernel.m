classdef J0Kernel < Kernel
    % Object of type kernel but optimized for J0 so that the computation of the 
    % radial quadrature goes faster
    % Radial quadrature is only made of 1 component !
    
    properties
        R
    end
    
    methods
        function[kernel] = J0Kernel(RR)
            kernel@Kernel(@(x)(besselj(0,RR*x)),@(x)(-RR*besselj(1,RR*x)))
            kernel = kernel.setScalFunc(@J0Kernel.scalFuncLaplace);
            kernel = kernel.setNormFunc(@J0Kernel.normFuncLaplace);
            kernel.PBounds = @(a,tol)([1, 1]);
            kernel.R = RR;
        end
        function[this] = dilatation(this,lambda)
            this = J0Kernel(lambda*this.R);
        end
    end
    
    
end

