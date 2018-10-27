function err = erreur(m, Uex, Uh, ef, type)
eps = 1e-6;
idx = [1:size(m.triangles,1)]';
% formule d'intégration dans Khat
[Xhat, omega] = formInt(2, 3);
xhat = Xhat(:,1); yhat = Xhat(:,2);
xA = m.vertices(m.triangles(idx,1),:);
xB = m.vertices(m.triangles(idx,2),:);
xC = m.vertices(m.triangles(idx,3),:);
x = kron(xA,(1 - xhat - yhat)) + kron(xB,xhat) + kron(xC,yhat);
w = kron(abs(m.detT(idx)),omega);
[n1,n2] = size(x);
eps1 = eps*ones(n1,1)*[1,0];
eps2 = eps*ones(n1,1)*[0,1];
Uapp = ef.ddl2int(m, idx, ef, Xhat)*Uh;
Uexact = Uex(x);

switch type
    case 'L2'
        err = sqrt(sum(w.*(Uexact-Uapp).^2));
    case 'H1'
        DxUapp = ef.ddl2int(m, idx, Dx(ef), Xhat)*Uh;
        DyUapp = ef.ddl2int(m, idx, Dy(ef), Xhat)*Uh;
        DxUexact = (Uex(x+eps1)-Uex(x-eps1))/(2*eps);
        DyUexact = (Uex(x+eps2)-Uex(x-eps2))/(2*eps);
        err = sqrt(sum(w.*((Uexact-Uapp).^2+(DxUexact-DxUapp).^2+(DyUexact-DyUapp).^2)));
    otherwise
        error('Unknown error type. Known types are ''L2'' or ''H1''.');
end
