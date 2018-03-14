function Dir = DirCond(ef, label, func)
DirDdl = find(~ismember(ef.lblDdl,label));
N2 = length(DirDdl);
N1 = ef.Nddl;
Dir.P = sparse(DirDdl,(1:N2)',ones(N2,1),N1,N2);
Dir.val = func(ef.ddl);
end
