function [ errH_12, errU1,errU0 ] = testN_OrderCV(curve,lambda, u0, Ns,varargin)
p = inputParser;
p.addOptional('quadNum',10);
p.addOptional('correcMethod','constantTerm');
p.addOptional('fullMatrix',false);
p.addOptional('a_factor',5);
p.addOptional('tol',1e-6);
p.addOptional('gmresTol',1e-10);
p.addOptional('prec',[]);
p.addOptional('omdxomu',[]);
p.KeepUnmatched = true;
p.parse(varargin{:});

errH_12 = zeros(length(Ns),1);
errU1 = zeros(length(Ns),1);
errU0 = zeros(length(Ns),1);

for i = 1:length(Ns)
    N = Ns(i);
    disp(N);
    meshAdapt = MeshCurve(curve,N,@cos,[-pi,0]);
    Vh =  weightedFEspace(meshAdapt,'P1','1/sqrt(1-t^2)',...
        'quadNum',p.Results.quadNum,'specialQuadSegs',1:meshAdapt.nseg);
    Wh = weightedFEspace(meshAdapt,'P1','sqrt(1-t^2)',...
        'quadNum',p.Results.quadNum,'specialQuadSegs',1:meshAdapt.nseg);
    omegaDxomega = Vh.omega_dx_omega;
    W = Vh.W;
    dM = omegaDxomega'*diag(W)*omegaDxomega;
    t1 = tic;
    Op_assemblingOptions = {'fullMatrix',p.Results.fullMatrix,'tol',p.Results.tol,'a_factor',p.Results.a_factor};
    Nw = hyperSingular_w(Vh,...
        'Op_opt',Op_assemblingOptions);
    Nwgalerk = Nw.galerkine;
    u0_h = Wh.secondMember(u0);
    if p.Results.fullMatrix
        lambda_h = variationalSol(Nwgalerk,u0_h,[],p.Results.gmresTol,N);
    else
        if isempty(p.Results.prec)
            lambda_h = variationalSol(Nwgalerk,u0_h,[],p.Results.gmresTol,N,Nwgalerk.concretePart);
        else
            lambda_h = variationalSol(Nwgalerk,u0_h,[],p.Results.gmresTol,N,prec);
        end
    end
    lambda_h = FE_func(Wh,lambda_h.fePart);
    t1 = toc(t1);
    fprintf('assembled and solved BEM equation in %s s \n',num2str(t1));
    errH_12(i) = sqrt((u0 | lambda - lambda_h));
    if ~isempty(p.Results.omdxomu)
        a1 = Vh.integral(p.Results.omdxomu.^2);
        a2 = sum(Vh.W.*(p.Results.omdxomu(Vh.gaussPoints)).*(omegaDxomega*lambda_h.v));
        a3 = lambda_h.v'*dM*lambda_h.v;
        errU1(i) = sqrt(a1 - 2*a2 + a3);
    end
    errU0(i) = sqrt((lambda - lambda_h | lambda - lambda_h));
end





end

