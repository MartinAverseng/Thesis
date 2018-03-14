classdef Y0Kernel < Kernel
    % Object of type kernel but optimized for Y0 so that the computation of the 
    % radial quadrature goes faster
    
    properties (Access = public)
        R; %such that kernel =  x -> Y0(Rx) 
    end
    
    methods
        function[kernel] = Y0Kernel(RR)
            func = @(x)(bessely(0,RR*x));
            der = @(x)(-RR*bessely(1,RR*x));
            kernel@Kernel(func,der);
            kernel.R = RR;
            kernel = kernel.setScalFunc(@(a,b,rho)(kernel.scalFuncHelmholtz(a,b,rho)));
            kernel = kernel.setNormFunc(@(a,b)(kernel.normFuncHelmholtz(a,b)));
            kernel = kernel.setStartFreq(RR);
        end
        
    end
    methods (Access = public)

        function[res] = scalFuncHelmholtz(this,a,b,rho)
            % Computes (a,b,\rho,R) -> \int_a^b x(Y_0(Rx))'J_1(\rho x) dx
            % ( = \int_a^b - Rx Y_1(Rx) J_1(\rho x) dx )
            RR = this.R;
            if ismember(RR,rho)
                error('bad choice of R ! It is not a root of Y0 or rho is not a root of J0')
            else
                res = helmholtzSP([a,b],rho,RR);
            end
        end
        function[res] = normFuncHelmholtz(this,a,b)
            % Computes (a,b) -> sqrt(2*pi*\int_a^b x [(Y0(Rx))']^2 dx))
            % (= sqrt(2*pi\int_a^b R^2 x Y_1(Rx)^2 dx)
            RR = this.R;
            T1 = b*RR*bessely(0,b*RR)^2;
            T2 = bessely(1,b*RR)*bessely(0,b*RR);
            T3 = b*RR*bessely(1,b*RR)^2;
            
            T4 = a*RR*bessely(0,a*RR)^2;
            T5 = bessely(1,a*RR)*bessely(0,a*RR);
            T6 = a*RR*bessely(1,a*RR)^2;
            
            res = sqrt(2*pi*RR/2*( b*(T1 - 2*T2 + T3) - a*(T4 - 2*T5 + T6)));
            
            % Explicit form
            %https://www.wolframalpha.com/input/?i=integral+from+a+to+b+of+R%5E2+x+*+(bessely(1,R+x))%5E2
        end
        function[rq] = radialQuadKernel(this,a,tol,varargin)
            % 1 Detect if R is a root of Y0, and otherwise, find next root
            if abs(bessely(0,this.R))>1e-12
                Rtemp = nextY0root(this.R);
            else
                Rtemp = this.R;
            end
            
            % Create a temporary kernel with the root
            Ktemp = Y0Kernel(Rtemp);
            % Dilatation factor
            lambda = this.R/Rtemp;
            atemp = a*lambda;
            rqTemp = RadialQuadrature(atemp,Ktemp,tol,varargin{:});
            % Use the previous radialQuad to form the final one : 
            rq = rqTemp.dilatation(lambda);            
        end
        
        
    end
    
    methods (Access = public)
        function[this] = dilatation(this,lambda)
            this = Y0Kernel(this.R*lambda);
        end
    end
end
