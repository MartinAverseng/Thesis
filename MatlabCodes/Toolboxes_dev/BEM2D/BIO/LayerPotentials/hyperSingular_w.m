classdef hyperSingular_w < BIO
    % weighted hypersingular operator N_omega
    properties
        r=1;
        k;
        ln_omega_reg;
        ln_omega2_reg
        correctionMethod = 'constantTerm';
    end
    
    methods
        function[this] = hyperSingular_w(VVh,varargin)
            assert(isa(VVh,'weightedFEspace'));
            assert(strcmp(VVh.weight_id,'1/sqrt(1-t^2)'));
            % We view the hypersingular opeerator as 
            % Nomega phi, psi = 
            % S1/omega omega dx omega phi , omega dx omega psi
            p = inputParser;
            p.addOptional('k',0);
            p.addOptional('Xdata',VVh.gaussPoints );
            p.addOptional('Op_opt',{})
            p.addOptional('r',1);
            p.parse(varargin{:});
            k = p.Results.k; XX = p.Results.Xdata; AopOpt = p.Results.Op_opt;
            r = p.Results.r;
            this.r = r;
            this.k = k;
            this.Vh = VVh;
            this.ln_omega_reg = this.Vh.ln_omega_reg(XX);
            
            if k > 0
                kern = 1i/4*H0Kernel(k);
                this.ln_omega2_reg = this.Vh.ln_omega2_reg(XX);
            else
                kern = (-1/(2*pi))*LogKernel(r);
                this.ln_omega2_reg = 0;
            end
            this.kernel = kern;
            this.Aop = Op(XX,kern,VVh.gaussPoints,AopOpt{:});
            this.AopOpt = AopOpt;
            this.X = XX;
            this.V = 'V';
        end
        function[] = set_X(this,XX,varargin)
%             if isequal(XX,this.X)
%                 return
%             end
            this.Aop = this.Aop.update_X(XX,varargin{:});
            this.ln_omega_reg = this.Vh.regularizeHypersingular(XX);
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
            
            M = Mkern*(AbstractMatrix.spdiag(this.Vh.W)*this.Vh.omega_dx_omega) + ...
                -1/(2*pi)*this.ln_omega_reg;
        end
        function[bili] = galerkine(this)
            Mkern = this.Aop;
            Wint = AbstractMatrix.spdiag(this.Vh.W).';
            M1 = Mkern*(Wint*this.Vh.omega_dx_omega) -1/(2*pi)*this.ln_omega_reg;
            mat1 = this.Vh.omega_dx_omega.'*Wint*M1;
            M2 = Mkern*(Wint*this.Vh.omega2) -1/(2*pi)*this.ln_omega2_reg;
            mat2 = this.Vh.omega2.'*Wint*M2;
            bili = BilinearForm(this.Vh,this.Vh,mat1 - this.k^2*mat2);
        end
    end
end


