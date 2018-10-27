N1 = 5e4;
N2 = 4e4;
rho = 3e-2;
dim = 3;
X = rand(N1,dim);
X(:,dim)=2*X(:,dim);
Y = rand(N2,dim);
Y(:,dim)=2*Y(:,dim)+0.5;
tic
B = myRangeSearch(rho,X,Y);
toc
tic
idx = rangesearch(X,Y,rho);
toc
disp('------VERIFICATION--------');
II = [];
JJ = [];
DD = [];
for i = 1:length(idx)
    A1 = idx{i};
    A2 = i*ones(size(A1));
    II = [II, A1];
    JJ = [JJ, A2];
    DD = [DD;sqrt(sum((X(A1',:)-Y(A2',:)).^2,2))];
end
D = sparse(II,JJ,DD',N1,N2);
if nnz(B-D)==0
    disp('OK')
else
    disp('ERROR')
end
%spy(B)