classdef HelmholtzPerturb < Kernel
    % Object of type kernel representing the real part of Gk - G0 where Gk is the kernel of
    % the Helmholtz equation with k frequency while G0 is the kernel of the
    % laplace operator. This kernel is analytical. 
    
    properties (Access = public)
        R; 
    end
    
    methods
        function[kernel] = HelmholtzPerturb(RR)
            kernel@Kernel(@func,@der);
            kernel.R = RR;
            kernel = kernel.setScalFunc(...
                @(a,b,rho)(...
                -1/4*kernel.scalFuncHelmholtz(a,b,rho))...
                +1/(2*pi)*kernel.scalFuncLog(a,b,rho)...
                );
            
            kernel = kernel.setNormFunc(@(a,b)(kernel.normFunc(a,b)));
            kernel = kernel.setStartFreq(RR);
            
            function[y] = func(x)
                y = 0*x;
                I = x>1e-10;
                J = x<=1e-10;
                y(I) = -1/4*(bessely(0,RR*x(I))) + 1/(2*pi)*log(RR*x(I));
                y(J) = 0.018451073777173;
            end
            
            function[y] = der(x)
                y = 0*x;
                I = x>1e-10;
                J = x<=1e-10;
                y(I) = 1/4*RR*bessely(1,RR*x(I)) + 1./(2*pi*x(I));
                y(J) = 0;
            end
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
        
        
        function[res] = normFunc(this,a,b)
            % Computes (a,b) -> \int_a^b x [-1/4.Y0(Rx))' + 1/(2.pi) ln(x)']^2 dx)
            RR = this.R;
            eugamma = double(eulergamma);
            prim = @(x)((1/4)*log(x)/pi^2+(1/32)*RR^2*x.^2.*(bessely(1, RR*x).^2-bessely(0, RR*x).*bessely(2, RR*x))-(1/4)*bessely(0, RR*x)/pi);
            prim0 = (1/8)*(1+2*log(2)-2*log(RR)-2*eugamma)/pi^2;
            
            if a~=0
                res = prim(b) - prim(a);
            else
                res = prim(b) - prim0;
            end
            
            % Explicit form
            %https://www.wolframalpha.com/input/?i=integral+from+a+to+b+of+R%5E2+x+*+(bessely(1,R+x))%5E2
        end
        
        
        
        function[res] = scalFuncLog(~,a,b,rho)
            % Computes (a,b,\rho) -> \int_a^b x(log(x))'J_1(\rho x) dx
            % ( = \int_a^b J_1(\rho x) dx )
            res = (besselj(0,rho*a) - besselj(0,rho*b))./rho;
            res = res(:)';
            % Explicit form :
            % https://www.wolframalpha.com/input/?i=integral+from+a+to+b+of++besselj(1,rho+x)
        end
        function[res] = normFuncLog(~,a,b)
            % Computes (a,b) -> \int_a^b x [(log(x))']^2 dx)
            % (= \int_a^b 1/x dx = log(b/a)
            res = sqrt(log(b/a));
        end
        
        
        
        
        
        
        function[rq] = radialQuadKernel(this,a,tol)
            % 1 Detect if R is a root of Y0, and otherwise, find next root
            if abs(bessely(0,this.R))>1e-12
                Rtemp = nextY0root(this.R);
            else
                Rtemp = this.R;
            end
            
            % Create a temporary kernel with the root
            Ktemp = HelmholtzPerturb(Rtemp);
            % Dilatation factor
            lambda = this.R/Rtemp;
            atemp = a*lambda;
            rqTemp = RadialQuadrature(atemp,Ktemp,tol);
            % Use the previous radialQuad to form the final one :
            rq = rqTemp.dilatation(lambda);
            
        end
        
        
    end
    
    methods (Access = public)
        function[this] = dilatation(this,lambda)
            this = HelmholtzPerturb(this.R*lambda);
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
