function out = iif(cond,a,b)
%IIF implements a ternary operator

% pre-assign out
out = b;

out(cond) = a(cond);