classdef Quad2D
    % Object that represents the approximation of a function in the form
    % f(x) \approx \sum_{\nu = 1}^{N_{\xi}} \hat{w}_{\nu}e^{i x \cdot \xi_{\nu}}
    % valid for a < |x| < 1, where \xi_{\nu} is a list of quadrature points
    % and \hat{w}_{\nu} a list of quadrature weights.
    properties (Constant,Hidden)
        % Constants developped to help choosing the number of points for
        % the circular discretisation of the spectrum. 
        gamma3 = 5;
        gamma4 = 0.33;
    end
    properties (Access = public)
        rq; % The radial quadrature
        w_nu; % Quadrature weights
        xi_nu; % Quadrature points
        time;
        offset;
    end
    
    methods
        function[q2d] = Quad2D(radialQuad)
            % Constructor
            ttime = tic;
            q2d.rq = radialQuad;
            [alpha,rho] = radialQuad.getAlpha;
            assert(rho(1)==0);
            
            q2d.offset = alpha(1);
            alpha = alpha(2:end);
            rho = rho(2:end);
            
            Ns = (fix((rho + Quad2D.gamma3 * rho.^(Quad2D.gamma4))/2)+1)*2;
            N = sum(Ns);
            xxi_nu = zeros(N,2);
            ww_nu = zeros(N,1);
            i1 = 1;
            % Loop on frequencies and add quadrature points as discretized
            % circle of radius this frequency. 
            for p = 1:length(rho)
                i2 = Ns(p) + i1 - 1;
                theta = (1:Ns(p))' * (2*pi)/Ns(p) ;
                xxi_nu(i1:i2,1) = rho(p)*cos(theta);
                xxi_nu(i1:i2,2) = rho(p)*sin(theta);
                ww_nu(i1:i2) = alpha(p)/Ns(p) * ones(Ns(p),1);
                i1 = i2 + 1;
            end
            q2d.w_nu = ww_nu;
            q2d.xi_nu = xxi_nu;
            q2d.time = toc(ttime);
        end
        function[] = disp(this)
            fprintf('2D quadrature of function %s \n',this.rq.getFunc)
            fprintf('Number of quadrature points : %s\n',num2str(length(this.w_nu)))
            [a,b] = this.rq.getInterval;
            fprintf('Domain of approximation : [%s, %s]\n',num2str(a),num2str(b))
            fprintf('Quadrature error : %s',this.rq.getQuadErrString);
        end
        function[el] = elem(this,X,Y,i,j)
            qj = zeros(size(X,1),1);
            qj(j) = 1;
            col = this.conv(X,Y,qj);
            el = col(i);
        end
        function[q,time] = conv(this,x,y,V)
            % This function computes \sum_{j=1}^N f(x(i)-y(j)) V(j)
            % where f is the function approximated by this object. 
            % This computes only accurately the far contribution, i.e. |x(i) - y(j)| > a
            % while the close interactions have to be corrected. 
            % |x(i) - y(j)| has to be <= 1 
            time = tic;
            xxi_nu = this.xi_nu;
            ww_nu = this.w_nu;
            Nxi = length(xxi_nu);
            tol = this.rq.getTol;
            %% Space to Fourier            
            V_nu = nufft2d3(size(y,1), y(:,1), y(:,2), ...
                V, -1, tol, Nxi, xxi_nu(:,1), xxi_nu(:,2) );
            %% Convolution becomes multiplication of Fourier weights
            fV_nu = V_nu.*ww_nu;            
            %% Back from Fourier to Space
            q = nufft2d3(Nxi, xxi_nu(:,1), xxi_nu(:,2), ...
                fV_nu, +1, tol, size(x,1), x(:,1), x(:,2));
            time = toc(time);
            %% Adding offset            
            q = q + this.offset*sum(V);
        end
        
    end
end

