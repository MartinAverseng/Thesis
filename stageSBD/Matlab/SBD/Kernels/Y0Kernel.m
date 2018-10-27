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
            kernel = kernel.setNormFunc(@(a,b,rho)(kernel.normFuncHelmholtz(a,b)));
            kernel = kernel.setPBounds(@Y0Kernel.PBounds);
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
            elseif a>0
                
                T1 = (a*RR)   *( bessely(0,a*RR) * besselj(1,a*rho) );
                T2 = (a*rho).*( bessely(1,a*RR) * besselj(0,a*rho) );
                T3 = (b*RR)   *( bessely(0,b*RR) * besselj(1,b*rho) );
                T4 = (b*rho).*( bessely(1,b*RR) * besselj(0,b*rho) );
                % allow rho to be a vector for vectorized computation !
                
                res = -RR./(RR^2 - rho.^2) .* (T1 - T2 - T3 + T4);
                % Explicit form :
                % https://www.wolframalpha.com/input/?i=integral+from+a+to+b+of++-R+x+bessely(1,R+x)+besselj(1,rho+x)
                res = res(:)';
            else
                T1 = pi*b*RR^2*besselj(1,b*rho)*bessely(0,b*RR);
                T2 = pi*b*RR*besselj(0,b*rho)*bessely(1,b*RR) + 2;
                % allow rho to be a vector for vectorized computation !
                
                res = 1./(pi*(RR^2-rho.^2)).*(T1 - rho.*T2);
                % Explicit form :
                % https://www.wolframalpha.com/input/?i=integral+from+0+to+b+of++-R+x+bessely(1,R+x)+besselj(1,rho+x)
                res = res(:)';
            end
        end
        function[res] = normFuncHelmholtz(this,a,b)
            % Computes (a,b) -> \int_a^b x [(Y0(Rx))']^2 dx)
            % (= \int_a^b R^2 x Y_1(Rx)^2 dx
            RR = this.R;
            T1 = b*RR*bessely(0,b*RR)^2;
            T2 = bessely(1,b*RR)*bessely(0,b*RR);
            T3 = b*RR*bessely(1,b*RR)^2;
            
            T4 = a*RR*bessely(0,a*RR)^2;
            T5 = bessely(1,a*RR)*bessely(0,a*RR);
            T6 = a*RR*bessely(1,a*RR)^2;
            
            res = sqrt(RR/2*( b*(T1 - 2*T2 + T3) - a*(T4 - 2*T5 + T6)));
            
            % Explicit form
            %https://www.wolframalpha.com/input/?i=integral+from+a+to+b+of+R%5E2+x+*+(bessely(1,R+x))%5E2
        end
        function[rq] = radialQuadKernel(this,a,~,tol)
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
            btemp = 1;
            rqTemp = RadialQuadrature(atemp,btemp,Ktemp,tol);
            % Use the previous radialQuad to form the final one : 
            rq = rqTemp.dilatation(lambda);
            
        end
        
        
    end
    
    methods (Access = public)
        function[this] = dilatation(this,lambda)
            this = Y0Kernel(this.R*lambda);
        end
    end
    
    methods (Static)
        function[loup] = PBounds(a,tol)
            % Only developped for the Laplace kernel for now... but often
            % good predictions
            loup = HeuristicBounds(a,tol);
            % If a warning was issued, do not be afraid, it's ok here. 
        end
    end
end
