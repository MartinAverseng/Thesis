classdef Somega < handle
    
    properties
        Vh@weightedFEspace; % The FE space on which the operator is defined
        k;
        Aop = Op;
        X;
        ln_reg;
    end
    
    methods
        function[this] = Somega(k,Vh,XX,varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            p.addOptional('r_var',1);
            p.parse(varargin{:});
            this.Vh = Vh;
            this.X = XX;
            this.ln_reg = this.Vh.regularize(X,'ln');
            if k > 0
                kern = 1i/4*H0Kernel(k);
            else
                kern = (-1/(2*pi))*LogKernel(p.Results.r_var);
            end
            this.Aop = Op(XX,kern,Vh.gaussPoints,varargin{:});
        end
        function[w] = weight(this)
            w = this.Vh.weight;
        end
        function[] = set_X(this,X,varargin)
            this.Aop = this.Aop.update_X(X,varargin{:});
            this.ln_reg = this.Vh.regularize(X,'ln');
        end
        function[M] = Mat(this)
            % This retunrs the matrix such that (M*U)_i = Su(x_i)
            Mkern = this.Aop;
            M = Mkern*(AbstractMatrix.spdiag(this.Vh.W)*this.Vh.phi) + ...
                -1/(2*pi)*this.ln_reg;
        end
        function[vals] = mtimes(S,lambda)
            assert(S.Vh.contains(lambda));
            vals = S.Mat*lambda.v;
        end
        function[reg] = special_reg(this)
            singInt = this.Vh.cell.singularIntegralDictionnary('ln');
            T = this.Vh.dofIndexes;
            [Ax,Bx,Ay,By] = this.Vh.mesh.edgesCoords;
            lnreg = regularize(this.Vh,this.X,'ln');
            segsToCorrec = this.Vh.I;
            dofsToCorrec = this.Vh.dofsOnSegments(segsToCorrec);
            % Proximity threshold, under which integrals are considered
            % singular.
            threshold = 2*max(this.Vh.mesh.length);
            I = this.Vh.dofCloseToPoints(this.X,threshold);
            dofsToLoop = intersect(find(~cellfun('isempty',I))',dofsToCorrec)';
            lines_correc = [];
            cols_correc = [];
            vals_correc = [];
            for dof_y = dofsToLoop
                lines_X = I{dof_y}';
                cols_dof = 0*lines_X + dof_y;
                X_k = this.X(lines_X,:);
                val_correc = zeros(length(lines_X),this.Vh.cell.Nb);
                for kk = 1:length(X_k)
                    x = X_k(kk,:);
                    for b = 1:this.Vh.cell.Nb;
                        seg = find(T(:,b)==dof_y);
                        assert(length(seg)<=1);
                        if isempty(seg)
                            val_correc(:,b) = 0;
                        else
                            A = [Ax(seg), Ay(seg)];
                            B = [Bx(seg), By(seg)];
                            u = (B-A)/norm(B-A);
                            [sA,sB] = this.Vh.mesh.s_seg(seg);
                            assert(abs((sB-sA)-norm(B-A))<1e-10);
                            y1 = @(s)(A(:,1) + (s-sA)*u(:,1));
                            y2 = @(s)(A(:,2) + (s-sA)*u(:,2));
                            r = @(s)(sqrt((x(:,1) - y1(s)).^2 + (x(:,2) - y2(s)).^2));
                            w = this.Vh.weight;
                            fun = @(s)(log(r(s)).*w(s));
                            [s_hat,w_hat] = gaussQuad(fun,1,sA,sB);
                            sref = (s_hat-sA)/(sB-sA);
                            phi_b = this.Vh.cell.basis_phi_b(b,sref);
                            val_correc(kk,b) = sum(w_hat.*phi_b)-singInt{b}(x,A,B);
                        end
                    end
                end
                lines_correc = [lines_correc;lines_X];%#ok
                cols_correc = [cols_correc;cols_dof];%#ok
                vals_correc = [vals_correc;sum(val_correc,2)];%#ok
            end
            reg = sparse(lines_correc,cols_correc,vals_correc,size(this.X,1),this.Vh.ndof);
            reg = lnreg + reg;
        end
    end
end


