function[K] = rigidite_simpleCouche(X,opt)

if nargin==1
    opt = 'P0';
end

if isequal(opt,'P0')
    p = 0;
elseif isequal(opt,'P1')
    p = 1;
else
    error('The value ''%s'' of arg ''opt'' is not supported',opt);
end

if p==0
    K = rigidite_simpleCoucheP0(X);
elseif p == 1
    K = rigidite_simpleCoucheP1(X);
else
    error('this error should never happen');
end
end


