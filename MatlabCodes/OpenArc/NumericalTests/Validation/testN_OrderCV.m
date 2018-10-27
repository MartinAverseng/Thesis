function [ errH_12, errL2 ] = testN_OrderCV(curve,lambda, u0, Ns,varargin)
p = inputParser;
p.addOptional('quadNum',10);
p.addOptional('correcMethod','constantTerm');
p.addOptional('fullMatrix',false);
p.addOptional('a_factor',5);
p.addOptional('tol',1e-6);
p.addOptional('gmresTol',1e-10);
p.addOptional('prec',[]);
p.KeepUnmatched = true;
p.parse(varargin{:});

errH_12 = zeros(length(Ns),1);
errL2 = zeros(length(Ns),1);

for i = 1:length(Ns)
    N = Ns(i);
    disp(N);
    meshAdapt = MeshCurve(curve,N,@cos,[-pi,0]);
    Vh =  weightedFEspace(meshAdapt,'P1','1/sqrt(1-t^2)',...
        'quadNum',p.Results.quadNum,'specialQuadSegs',1:meshAdapt.nseg);
    Wh = weightedFEspace(meshAdapt,'P1','sqrt(1-t^2)',...
        'quadNum',p.Results.quadNum,'specialQuadSegs',1:meshAdapt.nseg);
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
    errL2(i) = sqrt(real((lambda - lambda_h)|(lambda - lambda_h)));
end





end

