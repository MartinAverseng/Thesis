function u = solve(A, F, Dir)
if nargin == 2
    u=A\F;
elseif nargin == 3
    AD = Dir.P'*A*Dir.P;
    FD = Dir.P'*(F - A*Dir.val);
    uD = AD\FD;
    u = Dir.P*uD + Dir.val;
else
    error('Wrong number of arguments');
end
end
