classdef Op
    % Op objects represent a matrix A which elements are given by 
    % A(i,j) = G(X(i) - Y(j)) for some set of points X and Y. However, this 
    % object does not store te full matrix, instead it uses a Fourier 
    % approximation to reduce the amount of data stored and the time needed 
    % for computing a Matrix-Vector product. 
    properties
        rMax % when the user inputs X and Y as data set, rMax is an upper 
        % bound of the maximal distance between a point in X and a point in
        % Y. The data set are then rescaled by X/rMax and Y/rMax
        x,y; % rescaled data
        kernel% Kernel object. It contains the function G and other information needed for the method
        a,b % Bounds of the radial quadrature
        tol % required accuracy of the approximation 
        far % QuadD object, computes the far interaction 
        close % Sparse matrix for the close interactions + local corrections
        NCI,Nxi % Number of close interaction pairs and quadrature points 
        timesClose,timesFar,timeTotalAssembling % Timing info
    end
    
    methods (Access = public)
        %% Constructor
        function[op] = Op(XX,YY,kkernel,aa,bb,ttol,varargin)
            p = inputParser;
            p.addOptional('noCloseField',false) % Set to true to disable close field
            % This is useful when you know that clouds XX and YY are far
            % apart. 
            p.parse(varargin{:});
            vars = p.Results;
            noCloseField = vars.noCloseField;
            
            tTotalAssembling = tic;
            op.rMax = Op.rMaxCalc(XX,YY);
            
            % Rescaling clouds and kernel
            op.x = XX/op.rMax;
            op.y = YY/op.rMax;
            op.kernel = kkernel.dilatation(op.rMax);
            
            op.tol = ttol;
            op.a = aa;
            op.b = bb;
            
            % Create far and close opereators
            [ffar,NNxi,tFar] = op.kernel.quad2D(op.a,op.b,op.tol);
            [cclose,NNCI,tClose] = op.kernel.B(op.x,op.y,op.a,ffar,op.tol,noCloseField);
            
            op.far = ffar;
            op.close = cclose;
            op.Nxi = NNxi;
            op.NCI = NNCI;
            op.timesClose = tClose;
            op.timesFar = tFar;
            op.timeTotalAssembling = toc(tTotalAssembling);
        end
        
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
            
            % Matrix-vector time
            alinea =  '     ';
            
            q = randn(size(this,2),1);
            
            % Create a fake close matrix (in case this.close = all zero
            % sparse, for example when radiating to a far cloud of points)
            
            i = randsample(size(this,1),this.Nxi,true);
            j = randsample(size(this,2),this.Nxi,true);
            fakeClose = sparse(i,j,randn(this.Nxi,1),size(this,1),size(this,2));
            
            
            % Retrieve far and close operator
            ffar = this.far;
            cclose = this.close;
            
            % retrieve clouds
            xx = this.x;
            yy = this.y;
            
            % Chronoes : 
            % 1) far operator
            tFar = tic;
            C1 = ffar.conv(xx,yy,q);%#ok : Just for chrono
            tFar = toc(tFar);
            
            % 2) fake close operator
            tCloseFake = tic;
            C2 = fakeClose*q; %#ok : just for chrono
            tCloseFake = toc(tCloseFake);
            
            % 3) true close operator
            tCloseTrue = tic;
            C3 = cclose*q; %#ok : Just for chrono
            tCloseTrue = toc(tCloseTrue);
            
            % total time
            tTot = tCloseTrue + tFar;
            fprintf('- Matrix-vector product time : %.2f s \n',tTot)
            fprintf([alinea '- Close product : %.2f %% \n'],tCloseTrue/tTot*100)
            fprintf([alinea '- Far operator : %.2f %%\n'],tFar/tTot*100)
            if (tCloseTrue > tFar)
                str_for_a = 'reduce';
            else
                str_for_a = 'increase';
            end
            if or(tCloseTrue > 10*tFar,tCloseTrue < 1/10*tFar)
                fprintf('(for more balanced matrix-vector product time, %s parameter ''a'') \n',str_for_a);
            end
            
            % Estimate the time it would need to do the full product
            sparsityConstFakeB = length(i)/(size(this,1)*size(this,2));
            fullTimeGuess = tCloseFake/sparsityConstFakeB;
            
            propSBDvsFull = (tFar+tCloseTrue)/fullTimeGuess;
            fprintf('MV product approximated in about %.2f %% of full product time (estimated at %.2f s) \n',100*propSBDvsFull,fullTimeGuess);
            fprintf('Speed-up of %.2fx\n',1/propSBDvsFull);
            
        end
        function[str] = memorySize(this)
            N = this.Nxi + this.NCI;
            NumBytes = 8*N;
            str = Bytes2str(NumBytes);
        end
        
        %% mtimes function
        function[C] = mtimes(A,B)
            % Check dimensions
            assert(size(A,2) == size(B,1),'Impossible to perform matrix product, incompatible dimensions');
            % Pre-allocate
            C = zeros(size(A,1),size(B,2));
            ffar = A.far;
            cclose = A.close;
            xx = A.x;
            yy = A.y;
            for i = 1:size(B,2)
                q = B(:,i);
                C1= ffar.conv(xx,yy,q);
                C2 = cclose*q;
                C(:,i) = C1 + C2;
            end
            
        end
        
        %% Others
        
        function[el] = elem(A,i,j)
            qj = zeros(size(A,2),1);
            qj(j) = 1;
            col = A*qj;
            el = col(i);
            
        end
        function[vec] = size(this,i)
            if nargin ==1
                vec = [size(this.x,1),size(this.y,1)];
            elseif nargin ==2
                vec = [size(this.x,1),size(this.y,1)];
                vec = vec(i);
            else
                error('Bad number of input arguments')
            end
        end
        function[err,info] = validate(this,V,q,tVal,tol)
            tic
            yy = this.y;
            fun = this.kernel.func;
            
            xx = this.x(:,1);
            rxy = sqrt((xx(1)-yy(:,1)).^2 +(xx(2)-yy(:,2)).^2);
            tic; ignoredValue = sum(fun(rxy).*V); t = toc;  %#ok just to count time
            NBEM = min(fix(tVal/t),size(this,1));
            if NBEM == 0
                err = NaN;
                warning('Could not validate any sample on such a short amount of time')
            else
                valSample = sort(randsample(size(this,1),NBEM));
                q2 = zeros(NBEM,1);
                for i = 1:length(valSample)
                    id = valSample(i);
                    q2(i) = 0;
                    xx = this.x(id,:);
                    rxy = sqrt((xx(1)-yy(:,1)).^2 +(xx(2)-yy(:,2)).^2);
                    q2(i) = sum(fun(rxy).*V);
                end
                
                err = norm(q(valSample)-q2,'inf');
                t = toc;
                if err > tol
                    warning('Found samples with greater error than required accuracy')
                    info.badSamples = valSample(abs(q(valSample)-q2)>tol);
                    info.maxErr = err;
                else
                    propVal = NBEM/size(this,1)*100;
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
    end
    
end

