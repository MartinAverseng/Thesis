classdef Kernel
    % Kernel
    
    properties (SetAccess = protected, GetAccess = public)
        func,der,PBounds,scalFunc,normFunc,startFreq
    end
    
    methods
        function[kernel] = Kernel(func,der)
            if nargin == 0
                func = @(x)(0*x);
                der = @(x)(0*x);
            end                
            kernel.func = func;
            kernel.der = der;
            fun1 = @(rho)(@(x)(x.*der(x)*besselj(1,rho(:)'*x)));
            fun2 = @(x)(x.*der(x).^2);
            kernel.scalFunc = @(a,b,rho)(integral(fun1(rho),a,b,'ArrayValued',true));
            kernel.normFunc = @(a,b)(sqrt(integral(fun2,a,b)));
            kernel.PBounds = @(a,tol)([1, Inf]);
            kernel.startFreq = 0;
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
        function[rq] = radialQuadKernel(this,a,b,tol)
            rq = RadialQuadrature(a,b,this,tol);
        end
    end
    methods (Access = public)
        function[q2D,Nxi,times] = quad2D(this,a,b,tol)
            % Computes the 2D quadrature of the form 
            % \sum_{\nu = 1}^{N_{\xi}} \hat{w}_{\nu} e^{i x \cdot \xi_{\nu}}
            tTotalTime = tic;
            % First compute the radial quadrature
            rq = this.radialQuadKernel(a,b,tol);
            % Use it to compute the 2D quadrature
            q2D = Quad2D(rq);
            Nxi = length(q2D.w_nu);
            times.radialQuadrature = rq.times.total;
            times.quad2D = q2D.time;
            times.total = toc(tTotalTime);
        end
        function[BB,NCI,times] = B(this,x,y,a,q2D,tol,noCloseField)
            % Computes the close-range matrix associated to this kernel
            Nx = size(x,1);
            Ny = size(y,1);
            times.rangeSearch = 0;
            times.computeInteractions = 0;
            times.NUFFTcloseField = 0;
            times.assemble = 0;
            tTotalTime = tic;
            if or(a==0,noCloseField)
                NCI = 0;
                BB = sparse(Nx,Ny); % all zeros.
            else
                % Find close interactions
                tRangeSearch = tic;
                [I,rxy] = rangesearch(y,x,a);                
                jdx = cell2mat(I')';
                rxy = cell2mat(rxy')';
                idx = zeros(size(jdx));
                j = 1;
                for i=1:length(I);
                    idx(j:j+length(I{i})-1) = i;
                    j = j + length(I{i});
                end
                % Save memory
                clear I;
                times.rangeSearch = toc(tRangeSearch);
                NCI = length(rxy);
                Nxi = length(q2D.w_nu);
                
                % close-field
                if NCI ~= 0
                    tComputeInteractions = tic;
                    B1_inds = this.func(rxy);
                    B1_inds(or(isinf(B1_inds),isnan(B1_inds))) = 0;
                    
                    times.computeInteractions = toc(tComputeInteractions);
                    
                    % local correction of error due to NUFFT
                    tNUFFTcloseField = tic;
                    B2_inds = nufft2d3(Nxi, q2D.xi_nu(:,1), q2D.xi_nu(:,2), ...
                        q2D.w_nu, +1, tol, length(idx), x(idx,1) - y(jdx,1), x(idx,2) - y(jdx,2));
                    times.NUFFTcloseField = toc(tNUFFTcloseField);
                    
                    % Sparse matrix
                    tAssemble = tic;
                    BB_inds = B1_inds - B2_inds - q2D.offset;
                    BB = sparse(idx,jdx,BB_inds,Nx,Ny);
                    times.assemble = toc(tAssemble);
                else 
                    % No close interactions
                    BB = sparse(Nx,Ny); % all zeros
                end                
                
            end
            times.total = toc(tTotalTime);
            
        end
        
        function[this] = dilatation(old,lambda)
            oldFunc = old.func;
            oldDer = old.der;
            ffunc = @(x)(oldFunc(lambda*x));
            dder = @(x)(lambda*oldDer(lambda*x));
            this = Kernel(ffunc,dder);
            this = this.setStartFreq(lambda*old.startFreq);
            
        end
        
        
    end
end

