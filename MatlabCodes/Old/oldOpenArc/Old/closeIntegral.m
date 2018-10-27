function [ sp_mat] = closeIntegral(weightedVh,kernel,X,I)

TY = Vh.dofIndexes;
threshold = max(Vh.mesh.length);
I = Vh.dofCloseToPoints(X,threshold);
[YAx,YBx,YAy,YBy] = Vh.mesh.edgesCoords;

for dof_y = 1:nYdof
    lines_X = I{dof_y}';
    cols_dof = 0*lines_X + dof_y;
    X_k = X(lines_X,:);
    val_dof = zeros(length(lines_X),Nby);
    for k = 1:size(X_k,1)
        for b = 1:Nby;
            seg = find(TY(:,b)==dof_y);
            assert(length(seg)<=1);
            if isempty(seg)
                val_dof(:,b) = 0;
            else
                AY = [YAx(seg) YAy(seg)];
                BY = [YBx(seg) YBy(seg)];
                [sA,sB] = weightedVh.mesh.s_seg(seg);
                s = @(t)(sA + t*(sB - sA));
                u = (AY - BY);
                fun = @(t)(kernel.func(X - (AY + t*u)).*);
                val_dof(k,b) = singInt{b}(X_k,AY,BY);
            end
        end
    end
    
    
    lines_explicit = [lines_explicit;lines_X];%#ok
    cols_explicit = [cols_explicit;cols_dof];%#ok
    vals_explicit = [vals_explicit;sum(val_dof,2)];%#ok
    
end

explicit = sparse(lines_explicit,cols_explicit,vals_explicit,nX,nYdof);

sp_mat =  explicit - approx;

end





