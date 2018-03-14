classdef singleLayer < BIO
    
    properties
        r=1;
        k;
        ln_reg;
    end
    
    methods
        function[this] = singleLayer(k,VVh,XX,AopOpt,r)
            if ~exist('XX','var')||isempty(XX)
                XX = VVh.gaussPoints;
            end
            this.k = k;
            if ~exist('AopOpt','var')||isempty(AopOpt)
                AopOpt = {};
            end
            if ~exist('r','var')||isempty(r)
                r = 1;
            end
            this.r = r;
            this.Vh = VVh;
            this.ln_reg = this.Vh.regularize(XX,'ln');
            if k > 0
                kern = 1i/4*H0Kernel(k);
            else
                kern = (-1/(2*pi))*LogKernel(r);
            end
            this.kernel = kern;
            this.Aop = Op(XX,kern,VVh.gaussPoints,AopOpt{:});
            this.AopOpt = AopOpt;
            this.X = XX;
            this.V = 'V';
        end
        function[] = set_X(this,XX,varargin)
            if isequal(XX,this.X)
                return
            end
            this.Aop = this.Aop.update_X(XX,varargin{:});
            this.ln_reg = this.Vh.regularize(XX,'ln');
            this.X = XX;
        end
        function[this] = remesh(this,N)
            VVh = this.Vh.remesh(N);
            if isequal(this.X,this.Vh.gaussPoints)
                XX = VVh.gaussPoints;
            else
                XX = this.X;
            end
            this = singleLayer(this.k,VVh,XX,this.AopOpt,this.r);
        end
        function[M] = Mat(this)
            % This retunrs the matrix such that (M*U)_i = Su(x_i)
            Mkern = this.Aop;
            M = Mkern*(AbstractMatrix.spdiag(this.Vh.W)*this.testFunc) + ...
                -1/(2*pi)*this.ln_reg;
        end
    end
end


