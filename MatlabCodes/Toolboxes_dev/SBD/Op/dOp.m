classdef dOp
    % Op objects represent a matrix A which elements are given by
    % A(i,j) = grad(G)(X(i) - Y(j)) for some set of points X and Y. However, this
    % object does not store the full matrix, instead it uses a Fourier
    % approximation to reduce the amount of data stored and the time needed
    % for computing a Matrix-Vector product.
    properties
        dxOp@Op
        dyOp@Op
        kernel@Kernel = Kernel % Kernel object. It contains the function G and other information needed for the method
        tol % required accuracy of the approximation
        timesClose,timesFar,timeTotalAssembling % Timing info
        full = false;
    end
    
    methods (Access = public)
        %% Constructor
        function[dop] = dOp(X,k,Y,varargin)
            p = InputParser;
            p.addOptional('rMax',Op.rMaxCalc(X,Y))
            p.parse(varargin{:});
            optVars = p.Results;
            % Rescaling 
            
            Nx = size(X,1);
            Ny = size(Y,1);
            rMax = optVars.rMax;
            x = X/rMax;
            y = Y/rMax;
            k = k.dilatation(rMax);
            
            
            
            
            
            dop.dxOp = Op();
            N1 = size(X,1);
            N2 = size(Y,1);
            p = inputParser;
            p.KeepUnmatched = true;
            p.addOptional('noCloseField',false); % Set to true to disable close field
            % This is useful when you know that clouds XX and YY are far
            % apart.
            p.addOptional('noFarField',false);
            p.addOptional('fullMatrix',false);
            p.addOptional('a',Op.default_a(N1,N2));
            p.addOptional('a_factor',1);
            p.addOptional('tol',1e-3);
           
            p.parse(varargin{:});
            vars = p.Results;
            noFarField = vars.noFarField;
            noCloseField = vars.noCloseField;
            fullMatrix = vars.fullMatrix;
            aa = vars.a;
            a_factor = vars.a_factor;
            ttol = vars.tol;
            tTotalAssembling = tic;
            rMax = Op.rMaxCalc(X,Y)*1.1;
            % Rescaling clouds and kernel
            xx = X/rMax;
            yy = Y/rMax;
            k = k.dilatation(rMax);
            dop.kernel = k;
            if fullMatrix
                dop.a = 1;
                dop.tol = 0;
                noFarField = true;
                ffar1 = Quad2D; % empty Quad2D object.
                ffar2 = Quad2D; % empty Quad2D object.
                X1 = xx(:,1); X2 = xx(:,2);
                Y1 = yy(:,1); Y2 = yy(:,2);
                DIFF1 = repmat(X1,1,N2) - repmat(Y1',N1,1);
                DIFF2 = repmat(X2,1,N2) - repmat(Y2',N1,1);
                DIFFNorm = sqrt(DIFF1.^2 + DIFF2.^2);
                dop.concretePart = k.func(DIFFNorm);
                dop.concretePart(or(isnan(dop.concretePart),...
                    isinf(dop.concretePart))) = k.func(1e-8/dop.rMax);
                dop.far = ffar;
                dop.Nxi = NNxi;
                dop.NCI = N1*N2;
                dop.timesClose = toc(tTotalAssembling);
                dop.timesFar = 0;
                dop.timeTotalAssembling = dop.timesClose;
                dop.full = true;
            else
                dop.tol = ttol;
                dop.a = a_factor*aa;
                
                % Create far and close opereators
                [ffar,NNxi,tFar] = Op.quad2D(k,dop.a,dop.tol);
                if ~isa(k,'J0Kernel')
                    [cclose,NNCI,tClose] = Op.B(k,dop.x,dop.y,dop.a,ffar,dop.tol,noCloseField,dop.rMax);
                else
                    % No need to compute a close quadrature.
                    cclose = sparse([],[],[],N1,N2);
                    NNCI = 0;
                    tClose = 0;
                end
                if noFarField
                    [ffar,NNxi,tFar] = Op.quad2D(k,dop.a,dop.tol,noFarField); % empty one, to clean memory
                    % But for the close matrix, we need the step to compute the far
                    % quadrature.
                end
                dop.far = ffar;
                dop.concretePart = cclose;
                dop.abstractPart = @(vec)(ffar.conv(xx,yy,vec));
                dop.Nxi = NNxi;
                dop.NCI = NNCI;
                dop.timesClose = tClose;
                dop.timesFar = tFar;
                dop.timeTotalAssembling = toc(tTotalAssembling);
            end
            
            
        end
        %% Setters
        %% Display
        function[] = disp(this)
            
            printStarLine;
            fprintf('GENERAL INFO \n')
            fprintf('\nOperator of size %s x %s for kernel \n',num2str(size(this,1)),num2str(size(this,2)));
            disp(this.kernel.dilatation(1/this.rMax));
            
            rm = this.rMin;
            rM = this.rMax;
            fprintf('Close operator defined for r in [0, %.2f] \n',rm);
            fprintf('Far operator defined for [%.2f, %.2f] \n',rm,rM);
            fprintf('Number of close interactions : %s\n',num2str(this.NCI));
            fprintf('Number of radial quadrature points : %s\n',num2str(length(this.far.rq.rho)));
            fprintf('Number of quadrature points : %s\n',num2str(this.Nxi));
        end
        function[] = show(this)
            disp(this);
            showTime(this);
            showPerf(this);
            this.far.rq.show;
        end
        function[] = showTime(this)
            printStarLine;
            tClose = this.timesClose.total;
            tFar = this.timesFar.total;
            tTot = tClose + tFar;
            fprintf('OFFLINE PERFORMANCE \n')
            fprintf('\n- Operator assembling : %.2f s  :\n',tTot);
            alinea =  '     ';
            
            fprintf([alinea '- Close operator : %.2f s (%.2f %%): \n'],tClose,tClose/tTot*100);
            fprintf([alinea alinea '* Range search : %.2f %% \n'],this.timesClose.rangeSearch/tClose*100);
            fprintf([alinea alinea '* Computation of close interactions : %.2f %% \n'],this.timesClose.computeInteractions/tClose*100);
            fprintf([alinea alinea '* Local correction of the far-field (NUFFT) : %.2f %% \n'],this.timesClose.NUFFTcloseField/tClose*100);
            fprintf([alinea alinea '* Assembling close operator : %.2f %% \n'],this.timesClose.assemble/tClose*100);
            
            
            fprintf([alinea '- Far operator : %.2f s (%.2f %%): \n'],tFar,tFar/tTot*100);
            fprintf([alinea alinea '* Radial quadrature : %.2f %% \n'],this.timesFar.radialQuadrature/tFar*100);
            fprintf([alinea alinea '* 2D quadrature : %.2f %% \n'],this.timesFar.quad2D/tFar*100);
        end
        function[] = showPerf(this)
            printStarLine;
            fprintf('ONLINE PERFORMANCE \n')
            % Memory
            alinea =  '     ';
            fprintf('\n- Memory allocation : %s \n',this.memorySize)
            sparsityConst = this.NCI/(size(this,1)*size(this,2));
            fprintf([alinea '- Close operator : %.2f %% (%.2f %% sparsity)\n'],this.NCI/(this.Nxi + this.NCI)*100,100 - sparsityConst*100);
            fprintf([alinea '- Far operator : %.2f %%\n'],this.Nxi/(this.Nxi + this.NCI)*100)
            if (this.NCI > this.Nxi)
                str_for_a = 'reduce';
            else
                str_for_a = 'increase';
            end
            if or(this.NCI > 10*this.Nxi,this.NCI < 1/10*this.Nxi)
                fprintf('(for more balanced memory allocation, %s parameter ''a'') \n',str_for_a);
            end
            % Chrono of MV prod
            if (size(this,1)*size(this,2)<500000)
                tMV = 0;
                for i = 1:20
                    V = randn(this.N1,1);
                    tMVTic = tic;
                    q = this*V; %#ok
                    tMV = tMV + toc(tMVTic)/20;
                end
            elseif size(this,1)*size(this,2)<5000000
                tMV = 0;
                for i = 1:5
                    V = randn(this.N1,1);
                    tMVTic = tic;
                    q = this*V; %#ok
                    tMV = tMV + toc(tMVTic)/5;
                end
            else
                V = randn(this.N1,1);
                tMVTic = tic;
                q = this*V; %#ok
                tMV = toc(tMVTic);
            end
            fprintf([alinea '- MV product time : %s s \n'],num2str(tMV));
        end
        function[str] = memorySize(this)
            N = this.Nxi + this.NCI;
            NumBytes = 8*N;
            str = Bytes2str(NumBytes);
        end
        
        %% Others
        function[err,info] = validate(this,V,q,tVal,tol)
            if this.full
                error('nothing to validate');
            end
            totalTime = tic;
            yy = this.y;
            fun = this.kernel.func;
            
            xx = this.x(1,:);
            rxy = sqrt((xx(1)-yy(:,1)).^2 +(xx(2)-yy(:,2)).^2);
            tic; ignoredValue = sum(fun(rxy).*V); t = toc;  %#ok just to count time
            NBEM = min(fix(tVal/t/5),size(this,1));
            if NBEM == 0
                err = NaN;
                warning('Could not validate any sample on such a short amount of time')
            else
                shuffle = randperm(size(this,1));
                q2 = zeros(size(this,1),1);
                beginning = tic;
                try
                    for i = 1:size(this,1)
                        id = shuffle(i);
                        xx = this.x(id,:);
                        rxy = sqrt((xx(1)-yy(:,1)).^2 +(xx(2)-yy(:,2)).^2);
                        q2(i + max(toc(beginning)-tVal,0)*1i) = sum(fun(rxy).*V);
                    end
                catch
                    %Purposeful error 
                    i = i-1;
                end
                valSample = shuffle(1:i);
                err = norm(q(valSample)-q2(1:i),'inf');
                t = toc(totalTime);
                if err > tol
                    warning('Found samples with greater error than required accuracy')
                    info.badSamples = valSample(abs(q(valSample)-q2)>tol);
                    info.maxErr = err;
                else
                    propVal = min(i/size(this,1)*100,100);
                    fprintf('\n Validated %s%% of the samples in %s seconds \n',num2str(propVal),num2str(t));
                    info = [];
                end
                
                
            end
            
        end
        function[rm] = rMin(this)
            rm = this.a*this.rMax;
        end
        
    end
    
    methods (Static)
        function[rM] = rMaxCalc(X,Y)
            dist = @(x,x0) max( sqrt( (x(:,1)-x0(1)).^2 + (x(:,2)-x0(2)).^2) );
            rX   = dist(X,mean(X));
            rY   = dist(Y,mean(Y));
            rXY  = dist(mean(X),mean(Y));
            rM = rXY + rX  + rY;
        end
        function[aa] = default_a(N1,N2)
            aa = 1/(4*sqrt(sqrt(N1*N2)));
        end
        
        function[BB,NCI,times] = B(k,x,y,a,q2D,tol,noCloseField,rMax)
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
            else
                % Find close interactions
                tRangeSearch = tic;
                [I,rxy] = rangesearch(y,x,a*1.05);
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
                    B1_inds = k.func(rxy);
                    B1_inds(or(isinf(B1_inds),isnan(B1_inds))) = k.func(1e-8/rMax);
                    
                    times.computeInteractions = toc(tComputeInteractions);
                    
                    % local correction of error due to NUFFT
                    tNUFFTcloseField = tic;
                    B2_inds = nufft2d3(Nxi, q2D.xi_nu(:,1), q2D.xi_nu(:,2), ...
                        q2D.w_nu, +1, tol/10, length(idx), x(idx,1) - y(jdx,1), x(idx,2) - y(jdx,2));
                    times.NUFFTcloseField = toc(tNUFFTcloseField);
                    
                    % Sparse matrix
                    tAssemble = tic;
                    BB_inds = B1_inds - B2_inds - q2D.offset;
                    BB = sparse(idx(keepValues),jdx(keepValues),BB_inds(keepValues),Nx,Ny);
                    times.assemble = toc(tAssemble);
                else
                    % No close interactions
                    BB = sparse(Nx,Ny); % all zeros
                end
                
            end
            times.total = toc(tTotalTime);
            
        end
        function[q2D,Nxi,times] = quad2D(kernel,a,tol,noFarField)
            % Computes the 2D quadrature of the form
            % \sum_{\nu = 1}^{N_{\xi}} \hat{w}_{\nu} e^{i x \cdot \xi_{\nu}}
            tTotalTime = tic;
            if ~exist('noFarField','var')
                noFarField = false;
            end
            % First compute the radial quadrature
            if ~noFarField
                rq = kernel.radialQuadKernel(a,tol);
                % Use it to compute the 2D quadrature
                q2D = Quad2D(rq);
                Nxi = length(q2D.w_nu);
                times.radialQuadrature = rq.times.total;
                times.quad2D = q2D.time;
                times.total = toc(tTotalTime);
            else
                rq = RadialQuadrature;
                q2D = Quad2D(rq);
                Nxi = 0;
                times.radialQuadrature = 0;
                times.quad2D = 0;
                times.total = toc(tTotalTime);
            end
        end
    end
    
end

