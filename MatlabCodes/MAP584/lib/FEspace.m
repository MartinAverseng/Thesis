function ef = FEspace(dom, type, ordre)
%
switch type
    case 'Lagrange'
        ef = Lagrange(dom, ordre);
    otherwise
        error('FEspace.m : unknown type of finite element')
        
end
ef.op = 'Id';
end

