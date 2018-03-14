function [x,w,Xhat] = builtInteg(m, idx, quad)

if m.dim ==1
    % formule d'intégration dans Khat
    [Xhat, omega] = formInt(m.dim, quad);
    xA = m.vertices(m.elts(idx,1),:); xB = m.vertices(m.elts(idx,2),:);
    x = kron(xA,(1 - Xhat)) + kron(xB,Xhat);
    w = kron(m.vol_elts(idx),omega);
elseif m.dim == 2
    % formule d'intégration dans Khat
    [Xhat, omega] = formInt(m.dim, quad);
    xhat = Xhat(:,1); yhat = Xhat(:,2);
    xA = m.vertices(m.elts(idx,1),:);
    xB = m.vertices(m.elts(idx,2),:);
    xC = m.vertices(m.elts(idx,3),:);
    x = kron(xA,(1 - xhat - yhat)) + kron(xB,xhat) + kron(xC,yhat);
    w = kron(abs(m.vol_elts(idx)),omega);
end