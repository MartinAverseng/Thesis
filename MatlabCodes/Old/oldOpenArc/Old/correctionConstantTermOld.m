function[Mat] = correctionConstantTermOld(Vh,X,singKernel,I,varargin)
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
mesh = Vh.mesh;
ndof = Vh.ndof;
T = Vh.dofIndexes;
feCell = Vh.cell;
Nb = feCell.Nb;
gaussY = Vh.gaussPoints;
W = Vh.W;
func = singKernel.func;
[Ax,Bx,Ay,By] = mesh.edgesCoords; % edges [A, B] with A = [Ax, Ay], B = [Bx, By]

lines = []; cols = lines; vals = lines;
dofsToLoop = find(~cellfun('isempty',I))';
% Loop over Y dofs.
for i = dofsToLoop
    lines_X = I{i}';
    cols_dof = 0*lines_X + i;
    X_k = X(lines_X,:);
    M = size(X_k,1);
    correc_ki = zeros(length(lines_X),Nb);
    for b = 1:Nb;
        seg = find(T(:,b)==i);
        assert(length(seg)<=1);
        if isempty(seg)
            correc_ki(:,b) = 0;
        else
            A = repmat([Ax(seg) Ay(seg)],length(lines_X),1);
            B = repmat([Bx(seg) By(seg)],length(lines_X),1);
            cb = Vh.feCell.constantTerm(b,X,A,B);
            idx = (seg-1)*q + (1:q)';
            Yq = gaussY(idx,:);
            Wq = W(idx);
            
            rXkYq = sqrt((repmat(X_k(:,1),1,q) - repmat(Yq(:,1)',M,1)).^2 ...
                + (repmat(X_k(:,2),1,q) - repmat(Yq(:,2)',M,1)).^2);
            rXkYq(rXkYq < 1e-8) = 1e-8;
            Gkq = func(rXkYq);
            approxInt = Gkq*Wq;
            exactInt = singKernel.I0seg(X,AY,BY);
            correc_ki(:,b) = cb.*(exactInt - approxInt);
        end
    end
    
    lines = [lines;lines_X];%#ok
    cols = [cols;cols_dof];%#ok
    vals = [vals;sum(correc_ki,2)];%#ok
    
end
Mat = sparse(lines,cols,vals,M,ndof);
end