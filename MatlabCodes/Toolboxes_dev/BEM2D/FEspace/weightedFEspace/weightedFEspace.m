classdef weightedFEspace < FEspace
    %A space of functions defined on a mesh, endowed with a scalar product
    %of the form {f,g} = \int f(x)g(x)w(x) where w is a weight function.
    %A low order Gaussian quadrature for w is computed on each segment of
    %the mesh.
    
    properties
        I % indexes of the segments where to compute a special quadrature for w.
        weight % the weigth function defined in terms of the curvilinear abscissa
        weight_id = []; % If available
    end
    
    methods
        function[w] = FEweight(Vh,id)
            switch id
                case 'sqrt(1-t^2)'
                    L = sum(Vh.mesh.length);
                    w = @(s)(sqrt(s.*(L - s)));
                case '1/sqrt(1-t^2)'
                    L = sum(Vh.mesh.length);
                    w = @(s)(1./sqrt(s.*(L -s)));
                otherwise
                    error('unknown weight function');
                    
            end
        end
        function[this] = weightedFEspace(m,c,w,quadNum)
            this@FEspace(m,c,quadNum);
            q = this.quad.num;
            sVec = this.sVec;
            sVert = m.sVertices;
            if ischar(w)
                this.weight_id =w;
                this.weight = this.FEweight(w);
                w = this.weight;
            else
                this.weight = w;
                this.weight_id = [];
            end
            this.W = this.W.*w(sVec);
            T = this.dofIndexes;
            l = m.length;
            I = [1:5, m.nseg-(4:-1:0)];
            for i = 1:length(I)
                % quadrature
                seg_num = I(i);
                AB = m.segment(seg_num);
                A = AB(1); B = AB(2);
                s1 = sVert(A);
                s2 = sVert(B);
                A = this.mesh.vertex(A);
                B = this.mesh.vertex(B);
                
                [xhat,what] = gaussQuad(w,q,s1,s2);
                sref = (xhat - s1)/(s2-s1);
                idx = (seg_num-1)*q + (1:q)';
                u = (B-A);
                this.gaussPoints(idx,1) = A(1) + u(1)*sref;
                this.gaussPoints(idx,2) = A(2) + u(2)*sref;
                for b = 1:this.cell.Nb
                    [phi,dxphi] = this.cell.basis(sref);
                    dof = T(seg_num,b);
                    this.phi(idx,dof) = phi(:,b);
                    this.dphi(idx,dof) = dxphi(:,b)/l(seg_num);
                end
                this.W(idx) = what;
                this.sVec(idx) = s1 + sref*(s2 - s1);
            end
            Wdiag = AbstractMatrix.spdiag(this.W);
            M_mat = AbstractMatrix(this.phi'*Wdiag*this.phi);
            dM_mat = AbstractMatrix(this.dphi'*Wdiag*this.dphi);
            this.Mass = BilinearForm(this,this,M_mat);
            this.dMass = BilinearForm(this,this,dM_mat);
            this.selectedTest = this.phi;
            this.I = I;
            this.weight = w;
        end
        function[this] = remesh(this,N)
            m = this.mesh.remesh(N);
            if ~isempty(this.weight_id)
                this = weightedFEspace(m,this.cell.id,this.weight_id,this.quad.num);
            else
                this = weightedFEspace(m,this.cell.id,this.weight,this.quad.num);
            end
            
        end
        function[] = showWeight(this)
            s = this.sVec;
            plot(s,this.weight(s));
            xlabel('curvilinear absissa')
            title(sprintf('Weight function %s',this.weight_id));
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
            % for any list of M R   A2tildex = w0.*(Gkq*(Wq./this.weight(sq)));^2 points (X1,...XM) and for two R^2 points A and B
            
            M = size(X,1);
            me = this.mesh;
            L = sum(me.length);
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
                Ysing1 = [this.mesh.curve.x(-1) this.mesh.curve.y(-1)];
                Ysing2 = [this.mesh.curve.x(1) this.mesh.curve.y(1)];
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
                        [~,~,~,~,~,u,alpha,~,~] = parameters_singInt(X_k,repmat(A,size(X_k,1),1),repmat(B,size(X_k,1),1));
                        [sA,sB] = this.mesh.s_seg(seg);
                        s0 = sA -alpha;
                        s0(s0 < 0) = -s0(s0 < 0);
                        s0(s0>L) = 2*L-s0(s0 > L);
                        idx = (seg-1)*q + (1:q)';
                        Yq = gaussY(idx,:);
                        Wq = WY(idx);
                        w0 = this.weight(s0);
                        rXkYq = sqrt((repmat(X_k(:,1),1,q) - repmat(Yq(:,1)',Mk,1)).^2 ...
                            + (repmat(X_k(:,2),1,q) - repmat(Yq(:,2)',Mk,1)).^2);
                        rXkYq(rXkYq < 1e-8) = 1e-8;
                        Gkq = func(rXkYq);
                        sq = this.sVec(idx);
                        if ~ismember(seg,this.I)
                            % in this case, the only singularity comes
                            % from the kernel,
                            A2x = norm(B-A)*w0.*(singK.I0seg(X_k,A,B));
                            A2tildex = w0.*(Gkq*(Wq./this.weight(sq)));
                            correc_ki(:,b) = cb.*(A2x - A2tildex);
                        else
                            % In this case, we have to do things manually.
                            A2x1 = norm(B-A)*w0.*(singK.I0seg(X_k,A,B));
                            r1 = sqrt((X_k(:,1) - repmat(Ysing1(:,1)',Mk,1)).^2 ...
                            + (X_k(:,2) - repmat(Ysing1(:,2)',Mk,1)).^2);
                            r2 = sqrt((X_k(:,1) - repmat(Ysing2(:,1)',Mk,1)).^2 ...
                            + (X_k(:,2) - repmat(Ysing2(:,2)',Mk,1)).^2);
                            rsing = min(r1,r2);
                            A2x2 = norm(B-A)*w0.*func(rsing);
                            A2x = A2x1 - A2x2;
                            fsing = func(rsing);
                            A2tildex =(Gkq*Wq - fsing*sum(Wq));
                            for kk = 1:Mk
                                x = X_k(kk,:);
                                y1 = @(ss)(A(1) + (ss-sA)*u(kk,1));
                                y2 = @(ss)(A(2) + (ss-sA)*u(kk,2));
                                r = @(ss)(sqrt( ( x(1)-y1(ss) ).^2 + (x(2)-y2(ss)).^2));
                                fun = @(s)((func(r(s))-fsing(kk)).*(this.weight(s)-w0(kk)));
                                A2x(kk) = A2x(kk) + integral(fun,sA,sB,'AbsTol',1e-10);
                                
                            end
                            correc_ki(:,b) = cb.*(A2x - A2tildex);
                        end
                    end
                    
                end
                
                lines = [lines;lines_X];%#ok
                cols = [cols;cols_dof];%#ok
                vals = [vals;sum(correc_ki,2)];%#ok
                
            end
            Mat = sparse(lines,cols,vals,M,ndof);
        end
    end
end


