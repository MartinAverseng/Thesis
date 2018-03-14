classdef FEspace < handle
    % A FE space is an object defining an ensemble of functions defined on
    % a mesh and of the form f(x) = \sum \alpha_i \phi_i(x), where the
    % phi_i are basis functions described by the fe 'cell'.
    
    
    properties
        mesh@MeshCurve;
        cell@FEcell;
        quad@NumericalQuadrature
        phi
        dphi
        sVec
        gaussPoints
        W
        Mass@BilinearForm
        dMass@BilinearForm
        selectedTest
    end
    
    methods
        % Constructor
        function[this] = FEspace(m,c,quad)
            assert(ischar(c));
            switch c
                case 'P0'
                    this.cell = P0;
                case 'P1'
                    this.cell = P1;
                case 'P1disc'
                    this.cell = P1disc;
                case 'P2'
                    this.cell = P2;
                otherwise
                    error('unknown type of finite element space %s, or not supported yet',c);
            end
            this.mesh = m;
            if ~exist('quad','var')||isempty(quad)
                this.quad = inputOrDefaultQuad();
            else
                this.quad = inputOrDefaultQuad(quad);
            end
            
            [U,dU,ssVec,gaussPoints] = dof2point(this.mesh,this.cell,this.quad.xhat);
            this.phi = U;
            this.dphi = dU;
            this.sVec = ssVec;
            this.gaussPoints = gaussPoints;
            this.W = this.quad.weights(this.mesh);
            Wdiag = AbstractMatrix.spdiag(this.W);
            M_mat = AbstractMatrix(U'*Wdiag*U);
            dM_mat = AbstractMatrix(dU'*Wdiag*dU);
            this.Mass = BilinearForm(this,this,M_mat);
            this.dMass = BilinearForm(this,this,dM_mat);
            this.selectedTest = this.phi;
        end
        function[this] = remesh(this,N)
            m = this.mesh.remesh(N);
            this = FEspace(m,this.cell.id,this.quad.num);
            
        end
        function[N] = ndof(this)
            % number of degrees of freedom
            N = size(this.phi,2);
        end
        function[N] = nint(this)
            N = size(this.phi,1);
        end
        function[N] = nseg(this)
            % number of segments in the mesh
            N = this.mesh.nseg;
        end
        function[N] = nvert(this)
            % number of vertices in the mesh
            N = this.mesh.nvert;
        end
        function[T] = dofIndexes(this)
            T = this.cell.dofIndexes(this.mesh);
        end
        function[Z] = dofCoords(this)
            Z = this.cell.dofCoords(this.mesh);
        end
        function[fe_func] = Pi_h(this,func)
            % Creates the vector ndof x 1 of the coordinates of func on the
            % finite element space.
            % func must accept N x 2 arguments and return a N x 1 array
            Z = this.dofCoords;
            assert(isequal(size(Z),[this.ndof,2]));
            assert(isa(func,'R2toRfunc'),sprintf(...
                'Second argument should be of type ''R2toRfunc'' but received %s instead '...
                ,class(func)));
            v = func(Z);
            fe_func = FE_func(this,v);
        end
        function[out] = integral(this,u)
            u = func(u);
            if isempty(u.feSpace)
                u.feSpace = this;
            else
                assert(isequal(this,u.feSpace));
            end
            Wdiag = AbstractMatrix.spdiag(this.W);
            Phi = this.phi;
            gP = this.gaussPoints;
            [uh,f] = parts(u);
            out = sum(Wdiag*Phi*uh.v + Wdiag*f(gP));
        end
        function[l] = secondMember(this,k_func,varargin)
            k_func = func(this,k_func);
            [uh,f] = parts(k_func);
            vals = (f(this.gaussPoints).'.*this.W.')*this.phi;
            M = this.Mass;
            l = LinearForm(this,vals) + M*uh;
        end
        function[a] = meansBasis(this)
            a = sum(spdiags(this.W,0,length(this.W),length(this.W))*this.phi,1)';
        end
        function[dofs] = dofsOnSegments(this,segs)
            T = this.dofIndexes;
            dofs = T(segs,:);
            dofs = unique(dofs(:));
        end
        function[I] = dofCloseToPoints(this,X,threshold)
            I = isClose(X,this.dofCoords,threshold);
        end
        function[Mat] = exactCloseIntegrals(this,X,singKernel_id,I,varargin)
            % This function computes the matrix
            % Aij = \int_[this]K(x_i,y)phi_j(y) if I{j} contains i, 0
            % otherwise.
            % phi_j is the j-th basis function of the FE space or its
            % derivative if the argument 'derivative' is set to true.
            % K is the singular kernel represented by singKernel_id.
            p = inputParser;
            p.addOptional('derivative',false);
            p.parse(varargin{:});
            vars = p.Results;
            
            nX = size(X,1);
            m = this.mesh;
            ndof = this.ndof;
            T = this.dofIndexes;
            feCell = this.cell;
            Nb = feCell.Nb;
            
            [Ax,Bx,Ay,By] = m.edgesCoords; % edges [A, B] with A = [Ax, Ay], B = [Bx, By]
            if vars.derivative
                singInt = this.cell.singularIntegralDictionnaryDer(singKernel_id);
            else
                singInt = this.cell.singularIntegralDictionnary(singKernel_id);
            end
            
            lines = []; cols = lines; vals = lines;
            dofsToLoop = find(~cellfun('isempty',I))';
            % Loop over Y dofs.
            for dof_y = dofsToLoop
                lines_X = I{dof_y}';
                cols_dof = 0*lines_X + dof_y;
                X_k = X(lines_X,:);
                val_dofX = zeros(length(lines_X),Nb);
                for b = 1:Nb;
                    seg = find(T(:,b)==dof_y);
                    assert(length(seg)<=1);
                    if isempty(seg)
                        val_dofX(:,b) = 0;
                    else
                        AY = repmat([Ax(seg) Ay(seg)],length(lines_X),1);
                        BY = repmat([Bx(seg) By(seg)],length(lines_X),1);
                        val_dofX(:,b) = singInt{b}(X_k,AY,BY);
                    end
                end
                
                lines = [lines;lines_X];%#ok
                cols = [cols;cols_dof];%#ok
                vals = [vals;sum(val_dofX,2)];%#ok
                
            end
            Mat = sparse(lines,cols,vals,nX,ndof);
        end
        function[Mat] = approxCloseIntegrals(this,X,singKernel_id,I,varargin)
            % This function computes the matrix
            % Aij = \int_[this]K(x_i,y)phi_j(y) if I{j} contains i, 0
            % otherwise, with an approximate integral using the Gauss
            % points of this FE space.
            % phi_j is the j-th basis function of the FE space or its
            % derivative if the argument 'derivative' is set to true.
            % K is the singular kernel represented by singKernel_id.
            p = inputParser;
            p.addOptional('derivative',false);
            p.addOptional('constantTerm',false); % Set to true if you need
            % to correct only the constant term.
            p.parse(varargin{:});
            vars = p.Results;
            k = singularityDictionary(singKernel_id); % a kernel object.
            YGauss = this.gaussPoints;
            U = this.phi;
            
            WW = AbstractMatrix.spdiag(this.W);
            nX = size(X,1);
            ndof = this.ndof;
            if vars.derivative
                func = k.der;
            else
                func = k.func;
            end
            lines = []; cols = lines; vals = lines;
            dofsToLoop = find(~cellfun('isempty',I))';
            % Loop over Y dofs.
            for dof_y = dofsToLoop
                lines_X = I{dof_y}';
                cols_dof = 0*lines_X + dof_y;
                X_k = X(lines_X,:);
                l = dof_y;
                phiY_l = U(:,l);
                p_ind = phiY_l~=0;
                WY_p = WW.concretePart(p_ind,p_ind);
                if vars.constantTerm
                    UY_p = this.cell.constantTerm(X_k,dof_y); % This is a matrix.
                else
                    
                    Y_p_Gauss = YGauss(p_ind,:);
                    UY_p = phiY_l(p_ind);
                end
                n_p = nnz(phiY_l);
                n_k = length(lines_X);
                rXkYpGauss = sqrt((repmat(X_k(:,1),1,n_p) - repmat(Y_p_Gauss(:,1)',n_k,1)).^2 ...
                    + (repmat(X_k(:,2),1,n_p) - repmat(Y_p_Gauss(:,2)',n_k,1)).^2);
                rXkYpGauss(rXkYpGauss < 1e-8) = 1e-8;
                Gkp = func(rXkYpGauss);
                val_dofX = Gkp*(WY_p*UY_p);
                
                % Assembling matrix
                lines = [lines;lines_X];%#ok
                cols = [cols;cols_dof];%#ok
                vals = [vals;val_dofX];%#ok
                
            end
            Mat = sparse(lines,cols,vals,nX,ndof);
        end
        function[Mat] = regularize(this,X,singKernel_id,varargin)
            
            % Creates the Matrix  Mat such that Mat*U = \int_{\Gamma} kernel(\abs{X-y})u(y)dy -
            % \intApprox{\Gamma} kernel(\abs{X-y})u(y)dy
            % where \intApprox would be what one obtains with a simple numerical
            % quadrature, and u = sum_{i}U_i phi_i(y)
            % (this matrix is the corrective part accounting for singularities).
            
            p = inputParser;
            p.addOptional('threshold',nan);
            p.addOptional('derivative',false);
            p.addOptional('method','constantTerm');
            p.parse(varargin{:});
            vars = p.Results;
            threshold = vars.threshold;
            
            m = this.mesh;
            
            % searching where the integrals are singular
            if isnan(threshold)
                % Proximity threshold, under which integrals are considered
                % singular.
                threshold = 1.5*max(m.length);
            end
            
            I = this.dofCloseToPoints(X,threshold);
            
            if isequal(vars.method,'explicit')
                approx = this.approxCloseIntegrals(X,singKernel_id,I,'derivative',vars.derivative);
                explicit = this.exactCloseIntegrals(X,singKernel_id,I,'derivative',vars.derivative);
                Mat = explicit - approx;
            elseif isequal(vars.method,'constantTerm')
                Mat = this.correctionConstantTerm(X,singKernel_id,I,varargin);
            end
        end
        function[Mat] = correctionConstantTerm(this,X,singKernel_id,I,varargin)
            % Vh a weighted FE space.
            % This function computes the matrix
            % Aij = \int_[this]K(x_i,y)phi_j(y) if I{j} contains i, 0
            % otherwise.
            % phi_j is the j-th basis function of the FE space or its
            % derivative if the argument 'derivative' is set to true.
            % singKernel is the singular kernel. It must be able to return the handle
            % function corresponding to the singular kernel, and moreover must also
            % provide a function I(X,A,B) that computes \int_{[A,B]}\ln(|X - Y|)dY
            % for any list of M R^2 points (X1,...XM) and for two R^2 points A and B
            
            M = size(X,1);
            me = this.mesh;
            ndof = this.ndof;
            T = this.dofIndexes;
            feCell = this.cell;
            Nb = feCell.Nb;
            q = this.quad.num;
            gaussY = this.gaussPoints;
            WY = this.W;
            singK = singKernel(singKernel_id);
            func = singK.k.func;
            [Ax,Bx,Ay,By] = me.edgesCoords; % edges [A, B] with A = [Ax, Ay], B = [Bx, By]
            
            lines = []; cols = lines; vals = lines;
            dofsToLoop = find(~cellfun('isempty',I))';
            % Loop over Y dofs.
            for i = dofsToLoop
                lines_X = I{i}';
                cols_dof = 0*lines_X + i;
                X_k = X(lines_X,:);
                Mk = size(X_k,1);
                correc_ki = zeros(length(lines_X),Nb);
                for b = 1:Nb;
                    seg = find(T(:,b)==i);
                    assert(length(seg)<=1);
                    if isempty(seg)
                        correc_ki(:,b) = 0;
                    else
                        A = [Ax(seg) Ay(seg)];
                        B = [Bx(seg) By(seg)];
                        cb = this.cell.constantTerm(b,X_k,A,B);
                        idx = (seg-1)*q + (1:q)';
                        Yq = gaussY(idx,:);
                        Wq = WY(idx);
                        
                        rXkYq = sqrt((repmat(X_k(:,1),1,q) - repmat(Yq(:,1)',Mk,1)).^2 ...
                            + (repmat(X_k(:,2),1,q) - repmat(Yq(:,2)',Mk,1)).^2);
                        rXkYq(rXkYq < 1e-8) = 1e-8;
                        Gkq = func(rXkYq);
                        approxInt = Gkq*Wq;
                        exactInt = norm(B-A)*singK.I0seg(X_k,A,B);
                        correc_ki(:,b) = cb.*(exactInt - approxInt);
                    end
                end
                
                lines = [lines;lines_X];%#ok
                cols = [cols;cols_dof];%#ok
                vals = [vals;sum(correc_ki,2)];%#ok
                
            end
            Mat = sparse(lines,cols,vals,M,ndof);
        end
        function[I] = segmentsCloseToDofs(this,threshold)
            % Maybe to be removed
            hmax = max(this.mesh.length);
            Middles = this.mesh.middle;
            I = isClose(Middles,this.dofCoords,threshold + 1.5*hmax);
        end
        function[I] = TheirSegmentsCloseToMyDofs(me,them,threshold)
            % Maybe to be removed
            mymesh = me.mesh;
            myDofs = me.dofCoords;
            theirmesh = them.mesh;
            hmax = max(mymesh.length);
            theirMiddles = theirmesh.middles;
            I = isClose(theirMiddles,myDofs,threshold + 1.5*hmax);
        end
        function[I] = segmentsCloseToDofsOf(me,them,threshold)
            % Maybe to be removed
            mymesh = me.mesh;
            myMiddles = mymesh.middle;
            theirmesh = them.mesh;
            hmax = max(theirmesh.length);
            theirDofs = them.dofCoords;
            I = isClose(myMiddles,theirDofs,threshold + 1.5*hmax);
        end
        
        % Deprecated, to be removed
%         function[out] = BIO(this,varargin)
%             %% Input arguments parsing
%             p = inputParser;
%             p.addOptional('k','none');
%             p.addOptional('dk','none');
%             p.addOptional('varChange',@(Z)(Z)) %
%             p.addOptional('singularIntegral',NaN);
%             p.addOptional('singularIntegralDer',NaN);
%             p.addOptional('thresholdSingInt',NaN);
%             p.addOptional('refineQuad',20);
%             p.addOptional('compressSBD',true);
%             p.parse(varargin{:});
%             vars = p.Results;
%             varChange = vars.varChange;
%             singInt = vars.singularIntegral;
%             singIntDer = vars.singularIntegralDer;
%             threshold = vars.thresholdSingInt;
%             compressSBD = vars.compressSBD;
%             k = vars.k;
%             dk = vars.dk;
%             q = this.quad;
%             XGauss = this.gaussPoints;
%             quadRefine = NumericalQuadrature(vars.refineQuad);
%             numRefine = quadRefine.num;
%             if isequal(quadRefine.num,q.num)
%                 Urefine = this.U;
%                 dUrefine = this.dU;
%                 Wrefine = this.W;
%                 Xrefine = XGauss;
%             else
%                 [Urefine,dUrefine,~,Xrefine] = ...
%                     dof2point(this.mesh,this.cell,quadRefine.xhat);
%                 Wrefine = quadRefine.weights(this.mesh);
%             end
%             %% General variables.
%             Vh = this;
%             m = Vh.mesh;
%             refCell = Vh.cell;
%             
%             WW = this.W;
%             nint = size(XGauss,1);
%             nint_refine = size(Xrefine,1);
%             ndof = Vh.ndof;
%             T = Vh.dofIndexes;
%             Nb = refCell.Nb;
%             [x1,x2,y1,y2] = m.edgesCoords;
%             UU = this.U;
%             dUU = this.dU;
%             num = q.num;
%             %% Creation of the Kernel matrix
%             
%             if compressSBD
%                 tol = 1e-9;
%                 a = 1/(2*sqrt(nint));
%                 b = 1;
%                 if isequal(k,'none')
%                     Aop = AbstractMatrix(@(x)(0*x),nint,nint);
%                 else
%                     Aop = Op(varChange(XGauss),varChange(XGauss),k,a,b,tol);
%                 end
%                 if isequal(dk,'none')
%                     dAop = AbstractMatrix(@(x)(0*x),nint,nint);
%                 else
%                     dAop = Op(varChange(XGauss),varChange(XGauss),dk,a,b,tol);
%                 end
%             else
%                 X1 = XGauss(:,1); X2 = XGauss(:,2);
%                 DIFF1 = repmat(X1,1,nint) - repmat(X1',nint,1);
%                 DIFF2 = repmat(X2,1,nint) - repmat(X2',nint,1);
%                 DIFFNorm = sqrt(DIFF1.^2 + DIFF2.^2);
%                 if isequal(k,'none')
%                     Aop = AbstractMatrix(@(x)(0*x),nint,nint);
%                 else
%                     Mat = k.func(DIFFNorm);
%                     Mat(or(isnan(Mat),isinf(Mat)))=0;
%                     Aop = AbstractMatrix(Mat);
%                     
%                 end
%                 if isequal(dk,'none')
%                     dAop = AbstractMatrix(@(x)(0*x),nint,nint);
%                 else
%                     Mat = dk.func(DIFFNorm);
%                     Mat(or(isnan(Mat),isinf(Mat)))=0;
%                     dAop = AbstractMatrix(Mat);
%                 end
%             end
%             
%             %% Correction of singular integrals
%             correc = sparse([],[],[],nint,ndof);
%             fineQuad = sparse([],[],[],nint_refine,ndof);
%             
%             correcDer = sparse([],[],[],nint,ndof);
%             fineQuadDer = sparse([],[],[],nint_refine,ndof);
%             % Correction for the non-derivative part (the k part)
%             if ~isnan2(singInt)
%                 % * searching where the integrals are singular
%                 if isnan(threshold)
%                     % Proximity threshold, under which inegral are considered
%                     % singular.
%                     threshold = max(m.length)/2;
%                 end
%                 I = Vh.segmentsCloseToDofs(threshold);
%                 
%                 
%                 lines = [];
%                 cols = lines;
%                 vals = lines;
%                 lines_refine = [];
%                 cols_refine = lines_refine;
%                 vals_refine = lines_refine;
%                 for dof = 1:ndof
%                     neighbourSegs = I{dof};
%                     A1 = repmat((1:num)',1,length(neighbourSegs));
%                     A2 = repmat((neighbourSegs-1)*num,num,1);
%                     A1refine = repmat((1:numRefine)',1,length(neighbourSegs));
%                     A2refine = repmat((neighbourSegs-1)*numRefine,numRefine,1);
%                     lines_dof = A1(:) + A2(:);
%                     cols_dof = 0*lines_dof + dof;
%                     lines_dof_refine = A1refine(:) + A2refine(:);
%                     cols_dof_refine = 0*lines_dof_refine + dof;
%                     X_k = XGauss(lines_dof,:);
%                     X_krefine = Xrefine(lines_dof_refine,:);
%                     vals_dof_refine = zeros(length(lines_dof_refine),Nb);
%                     
%                     for b = 1:Nb;
%                         seg = find(T(:,b)==dof);
%                         assert(length(seg)<=1);
%                         if isempty(seg)
%                             vals_dof_refine(:,b) = 0;
%                         else
%                             A = repmat([x1(seg) y1(seg)],length(lines_dof_refine),1);
%                             B = repmat([x2(seg) y2(seg)],length(lines_dof_refine),1);
%                             vals_dof_refine(:,b) = singInt{b}(X_krefine,A,B);
%                         end
%                     end
%                     
%                     % Old values that need to be substracted :
%                     l = dof;
%                     phi_l = UU(:,l);
%                     p_ind = phi_l~=0;
%                     X_p = XGauss(p_ind,:);
%                     W_p = WW(p_ind);
%                     U_p = phi_l(p_ind);
%                     n_p = nnz(phi_l);
%                     n_k = length(lines_dof);
%                     rXkXp = sqrt((repmat(X_k(:,1),1,n_p) - repmat(X_p(:,1)',n_k,1)).^2 ...
%                         + (repmat(X_k(:,2),1,n_p) - repmat(X_p(:,2)',n_k,1)).^2);
%                     Gkp = k.func(rXkXp);
%                     Gkp(or(isinf(Gkp),isnan(Gkp))) = 0;
%                     valCorrec = Gkp*(W_p.*U_p);
%                     
%                     % Assembling matrix
%                     lines = [lines;lines_dof]; % Is it possible to know the
%                     % size of lines, cols and vals in advance ?
%                     cols = [cols;cols_dof];
%                     vals = [vals;-valCorrec];
%                     
%                     lines_refine = [lines_refine;lines_dof_refine];
%                     cols_refine = [cols_refine;cols_dof_refine];
%                     vals_refine = [vals_refine;sum(vals_dof_refine,2)];
%                     
%                 end
%                 
%                 correc = sparse(lines,cols,vals,nint,ndof);
%                 fineQuad = sparse(lines_refine,cols_refine,vals_refine,nint_refine,ndof);
%             end
%             % Correction for the derivative part (the dk part)
%             if ~isnan2(singIntDer)
%                 % * searching where the integrals are singular
%                 if isnan(threshold)
%                     % Proximity threshold, under which inegral are considered
%                     % singular.
%                     threshold = max(m.length)/2;
%                 end
%                 I = Vh.segmentsCloseToDofs(threshold);
%                 
%                 
%                 lines = [];
%                 cols = lines;
%                 vals = lines;
%                 lines_refine = [];
%                 cols_refine = lines_refine;
%                 vals_refine = lines_refine;
%                 for dof = 1:ndof
%                     neighbourSegs = I{dof};
%                     A1 = repmat((1:num)',1,length(neighbourSegs));
%                     A2 = repmat((neighbourSegs-1)*num,num,1);
%                     A1refine = repmat((1:numRefine)',1,length(neighbourSegs));
%                     A2refine = repmat((neighbourSegs-1)*numRefine,numRefine,1);
%                     lines_dof = A1(:) + A2(:);
%                     cols_dof = 0*lines_dof + dof;
%                     lines_dof_refine = A1refine(:) + A2refine(:);
%                     cols_dof_refine = 0*lines_dof_refine + dof;
%                     X_k = XGauss(lines_dof,:);
%                     X_krefine = Xrefine(lines_dof_refine,:);
%                     vals_dof_refine = zeros(length(lines_dof_refine),Nb);
%                     
%                     
%                     for b = 1:Nb;
%                         seg = find(T(:,b)==dof);
%                         assert(length(seg)<=1);
%                         if isempty(seg)
%                             vals_dof_refine(:,b) = 0;
%                         else
%                             A = repmat([x1(seg) y1(seg)],length(lines_dof_refine),1);
%                             B = repmat([x2(seg) y2(seg)],length(lines_dof_refine),1);
%                             vals_dof_refine(:,b) = singIntDer{b}(X_krefine,A,B);
%                         end
%                     end
%                     
%                     % Old values that need to be substracted :
%                     l = dof;
%                     dphi_l = dUU(:,l);
%                     p_ind = dphi_l~=0;
%                     X_p = XGauss(p_ind,:);
%                     W_p = WW(p_ind);
%                     dU_p = dphi_l(p_ind);
%                     n_p = nnz(dphi_l);
%                     n_k = length(lines_dof);
%                     rXkXp = sqrt((repmat(X_k(:,1),1,n_p) - repmat(X_p(:,1)',n_k,1)).^2 ...
%                         + (repmat(X_k(:,2),1,n_p) - repmat(X_p(:,2)',n_k,1)).^2);
%                     Gkp = dk.func(rXkXp);
%                     Gkp(rXkXp<1e-10) = dk.func(1e-10);
%                     valCorrec = Gkp*(W_p.*dU_p);
%                     
%                     % Assembling matrix
%                     lines = [lines;lines_dof]; % Is it possible to know the
%                     % size of lines, cols and vals in advance ?
%                     cols = [cols;cols_dof];
%                     vals = [vals;-valCorrec];
%                     
%                     lines_refine = [lines_refine;lines_dof_refine];
%                     cols_refine = [cols_refine;cols_dof_refine];
%                     vals_refine = [vals_refine;sum(vals_dof_refine,2)];
%                     
%                 end
%                 correcDer = sparse(lines,cols,vals,nint,ndof);
%                 fineQuadDer = sparse(lines_refine,cols_refine,vals_refine,nint_refine,ndof);
%             end
%             %% Creation of the bilinear form
%             mv_p = @(x)(...
%                 real(...
%                 UU'*(WW.*(Aop.mv_prod(WW.*(UU*x)) + correc*x))...
%                 +	Urefine'*(Wrefine.*(fineQuad*x)))...
%                 + real(...
%                 dUU'*(WW.*(dAop.mv_prod(WW.*(dUU*x)) + correcDer*x))...
%                 + dUrefine'*(Wrefine.*(fineQuadDer*x)))...
%                 );
%             out = BilinearForm(this,mv_p);
%         end
        % Deprecated, to be removed
        function[S] = singleLayerPotential(this,k,varargin)
            
            p = inputParser;
            p.addOptional('r',1);
            p.addOptional('compressSBD',true);
            p.parse(varargin{:});
            vars = p.Results;
            r = vars.r;
            compressSBD = vars.compressSBD;
            if nargin == 1
                k=0;
            elseif k~=0
                if r~=1
                    warning('value of parameter ''r'' ignored because k~=0')
                end
                r = k;
            end
            assert(r>0);
            G0 = LogKernel(1);
            singInt = this.cell.singularIntegralDictionnary('ln');
            a = this.meansBasis();
            S0 = -1/(2*pi)*this.BIO('k',G0,'singularIntegral',singInt,'compressSBD',compressSBD) + -1/(2*pi)*log(r)*AbstractMatrix(@(x)(a*(a'*x)),this.ndof,this.ndof);
            if k ~= 0
                Gk_re = HelmholtzPerturb(k);
                Gk_im = J0Kernel(k);
                Sk_re = this.BIO('k',Gk_re,'compressSBD',compressSBD);
                Sk_im = 1/4*this.BIO('k',Gk_im,'compressSBD',compressSBD);
                S = S0 + Sk_re + 1i*Sk_im;
            else
                S = S0;
            end
        end
        % Deprecated, to be removed
        function[T] = doubleLayerPotential(this,varargin)
            % For now, I don't know how to compress the kernel, so full
            % version.
            p = inputParser;
            p.addOptional('refineQuad',6);
            p.addOptional('thresholdSingInt',NaN);
            p.parse(varargin{:});
            vars = p.Results;
            X = this.zVec;
            q = this.quad;
            threshold = vars.thresholdSingInt;
            quadRefine = NumericalQuadrature(vars.refineQuad);
            numRefine = quadRefine.num;
            if isequal(quadRefine.num,q.num)
                Urefine = this.U;
                Wrefine = this.W;
                Xrefine = X;
            else
                [Urefine,~,~,Xrefine] = ...
                    dof2point(this.mesh,this.cell,quadRefine.xhat);
                Wrefine = quadRefine.weights(this.mesh);
            end
            Vh = this;
            m = Vh.mesh;
            refCell = Vh.cell;
            WW = this.W;
            nint = size(X,1);
            nint_refine = size(Xrefine,1);
            ndof = Vh.ndof;
            Tab = Vh.dofIndexes;
            Nb = refCell.Nb;
            [x1,x2,y1,y2] = m.edgesCoords;
            UU = this.U;
            num = q.num;
            
            
            X1 = X(:,1);
            X2 = X(:,2);
            DIFF1 = repmat(X1,1,nint) - repmat(X1',nint,1);
            DIFF2 = repmat(X2,1,nint) - repmat(X2',nint,1);
            DIFFNorm = sqrt(DIFF1.^2 + DIFF2.^2);
            A1 = DIFF1./DIFFNorm.^2;
            A2 = DIFF2./DIFFNorm.^2;
            A1(or(isnan(A1),isinf(A1))) = 0;
            A2(or(isnan(A2),isinf(A2))) = 0;
            N = m.normalVector;
            N = N./repmat(sqrt(N(:,1).^2 + N(:,2).^2),1,2);
            N1 = spdiags(reshape(repmat(N(:,1)',num,1),nint,1),0,nint,nint);
            % Replicate the normal vector (first component) for each integration point.
            % And make a sparse diag out of it.
            N2 = spdiags(reshape(repmat(N(:,2)',num,1),nint,1),0,nint,nint);
            Mat = A1*N1+ A2*N2;
            Aop = AbstractMatrix(Mat);
            %% Correction of singular integrals
            
            singInt = refCell.singularIntegralDictionnary('doubleLayer');
            
            % Correction for the non-derivative part (the k part)
            
            if isnan(threshold)
                % Proximity threshold, under which inegral are considered
                % singular.
                threshold = max(m.length)/2;
            end
            I = Vh.segmentsCloseToDofs(threshold);
            
            lines = [];
            cols = lines;
            vals = lines;
            lines_refine = [];
            cols_refine = lines_refine;
            vals_refine = lines_refine;
            for dof = 1:ndof
                neighbourSegs = I{dof};
                A1 = repmat((1:num)',1,length(neighbourSegs));
                A2 = repmat((neighbourSegs-1)*num,num,1);
                A1refine = repmat((1:numRefine)',1,length(neighbourSegs));
                A2refine = repmat((neighbourSegs-1)*numRefine,numRefine,1);
                lines_dof = A1(:) + A2(:);
                cols_dof = 0*lines_dof + dof;
                lines_dof_refine = A1refine(:) + A2refine(:);
                cols_dof_refine = 0*lines_dof_refine + dof;
                X_krefine = Xrefine(lines_dof_refine,:);
                vals_dof_refine = zeros(length(lines_dof_refine),Nb);
                
                for b = 1:Nb;
                    seg = find(Tab(:,b)==dof);
                    assert(length(seg)<=1);
                    if isempty(seg)
                        vals_dof_refine(:,b) = 0;
                    else
                        A = repmat([x1(seg) y1(seg)],length(lines_dof_refine),1);
                        B = repmat([x2(seg) y2(seg)],length(lines_dof_refine),1);
                        Nloc = repmat(N(seg,:),length(lines_dof_refine),1);
                        vals_dof_refine(:,b) = singInt{b}(X_krefine,A, B,Nloc);
                    end
                end
                
                % Old values that need to be substracted :
                l = dof;
                phi_l = UU(:,l);
                p_ind = phi_l~=0;
                W_p = WW(p_ind);
                U_p = phi_l(p_ind);
                
                Gkp = Mat(lines_dof,p_ind);
                valCorrec = Gkp*(W_p.*U_p);
                
                % Assembling matrix
                lines = [lines;lines_dof]; % Is it possible to know the
                % size of lines, cols and vals in advance ?
                cols = [cols;cols_dof];
                vals = [vals;-valCorrec];
                
                lines_refine = [lines_refine;lines_dof_refine];
                cols_refine = [cols_refine;cols_dof_refine];
                vals_refine = [vals_refine;sum(vals_dof_refine,2)];
                
            end
            correc = sparse(lines,cols,vals,nint,ndof);
            fineQuad = sparse(lines_refine,cols_refine,vals_refine,nint_refine,ndof);
            
            %% Creating Bilinear form
            mv_p = @(x)(...
                real(...
                UU'*(WW.*(Aop.mv_prod(WW.*(UU*x)) + correc*x))...
                +	Urefine'*(Wrefine.*(fineQuad*x))));
            
            T = 1/pi*BilinearForm(Vh,mv_p);
            
        end
        % Deprecated, to be removed
        function[R] = hyperSingularOperator(this,varargin)
            p = inputParser;
            p.addOptional('zeroMeanModif',0)
            p.parse(varargin{:});
            vars = p.Results;
            z = vars.zeroMeanModif;
            
            G = LogKernel(1);
            singIntDer = this.cell.singularIntegralDictionnaryDer('ln');
            R = -1/(2*pi)*this.BIO('dk',G,'singularIntegralDer',singIntDer,'compressSBD',compressSBD);
            if z~= 0
                a = this.meansBasis();
                R = R + z*AbstractMatrix(@(x)(a*(a'*x)),this.ndof,this.ndof);
            end
        end
        
        function[bool] = contains(this,uh)
            if isa(uh,'FE_func')
                bool = isequal(uh.feSpace,this);
            else
                bool = false;
            end
        end
        
    end
end