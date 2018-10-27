function A = intg(dom, varargin)
global mesh
p = inputParser;
p.addParamValue('quad',3);
p.addParamValue('func',[]);
p.addParamValue('test',[]);
p.addParamValue('trial',[]);
p.addParamValue('label',-1);
p.parse(varargin{:});
res = p.Results;

if res.label ~= -1
    switch dom.dim
        case 0
            res.label = unique(mesh.lab_points);
        case 1
            res.label = unique(mesh.lab_edges);
        case 2
            res.label = unique(mesh.lab_triangles);
        case 3
            res.label = unique(mesh.lab_tetras);
        otherwise
            error('Incompatible dimension of domain');
    end
end
switch dom.dim
    case 0
        idx = find(ismember(mesh.lab_points, res.label));
        error('Not implemented');
    case 1
        idx = find(ismember(mesh.lab_edges, res.label));
        % formule d'intégration dans Khat
        [Xhat, omega] = formInt(dom.dim, res.quad);
        xA = mesh.vertices(mesh.edges(idx,1),:);
        xB = mesh.vertices(mesh.edges(idx,2),:);
        x = kron(xA,(1 - Xhat)) + kron(xB,Xhat);
        w = kron(mesh.vol_edges(idx),omega);
    case 2
        idx = find(ismember(mesh.lab_triangles, res.label));
        % formule d'intégration dans Khat
        [Xhat, omega] = formInt(dom.dim, res.quad);
        xhat = Xhat(:,1); yhat = Xhat(:,2);
        xA = mesh.vertices(mesh.triangles(idx,1),:);
        xB = mesh.vertices(mesh.triangles(idx,2),:);
        xC = mesh.vertices(mesh.triangles(idx,3),:);
        x = kron(xA,(1 - xhat - yhat)) + kron(xB,xhat) + kron(xC,yhat);
        w = kron(abs(mesh.vol_triangles(idx)),omega);
    case 3
        idx = find(ismember(mesh.lab_tetras, res.label));
        % formule d'intégration dans Khat
        [Xhat, omega] = formInt(dom.dim, res.quad);
        xhat = Xhat(:,1); yhat = Xhat(:,2); zhat = Xhat(:,3);
        xA = mesh.vertices(mesh.tetras(idx,1),:);
        xB = mesh.vertices(mesh.tetras(idx,2),:);
        xC = mesh.vertices(mesh.tetras(idx,3),:);
        xD = mesh.vertices(mesh.tetras(idx,4),:);
        x = kron(xA,(1 - xhat - yhat - zhat)) + kron(xB,xhat) + kron(xC,yhat) + kron(xD,zhat);
        w = kron(abs(mesh.vol_tetras(idx)),omega);
    otherwise
        error('Incompatible dimension of domain');
end

if isa(res.func, 'function_handle')
    w = w.*res.func(x);
end

if (~isempty(res.test))
    if (~isempty(res.trial))% Forme bilineaire
        Nint = size(w,1);
        W = spdiags(w,0,Nint,Nint);
        Atest = res.test.ddl2int(dom, idx, Xhat, res.test);
        Atrial = res.trial.ddl2int(dom, idx, Xhat, res.trial);
        A = Atest'*W*Atrial;
    else % Forme lineaire
        Atest = res.test.ddl2int(dom, idx, Xhat, res.test);
        A = Atest'*w;
    end
else
    A = sum(w);
end
end