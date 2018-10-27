function [ D ] = doubleLayer( k,Xh,Vh,varargin)
p = inputParser;
p.addOptional('r',1);
p.addOptional('fullMatrix',false);
p.addOptional('refineQuad',10)
p.addOptional('tol',1e-3)
p.addOptional('a_factor',1);
p.addOptional('noFarField',false);
p.parse(varargin{:});
vars = p.Results;
r = vars.r;
refineQuad = vars.refineQuad;
tol = vars.tol;

assert(and(isa(Xh,'FEspace'),isa(Yh,'FEspace')),'The second and third arguments must be FEspace Objects')
if k==0
    logK = LogKernel(r);
    S = BIOGalerkine(Xh,Yh,'U',logK,'V','tol',tol,varargin{:});
    S = -1/(2*pi)*(S + regularize(Xh,Yh,'U','ln','V','refineQuad',refineQuad));
else
    J0 = J0Kernel(k);
    Y0 = Y0Kernel(k);
    J0BEM = BIOGalerkine(Xh,Yh,'U',J0,'V','tol',tol,varargin{:});
    Y0BEM = BIOGalerkine(Xh,Yh,'U',Y0,'V','tol',tol,varargin{:});
    % Y0(z) ~ 2/pi*ln(z) + ...
    Y0BEM = Y0BEM + 2/pi*regularize(Xh,Yh,'U','ln','V','refineQuad',refineQuad);
    S = 1i/4*J0BEM - 1/4*Y0BEM;
end



end

