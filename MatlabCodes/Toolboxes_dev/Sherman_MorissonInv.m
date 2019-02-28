% Shermann-Morisson backslash
function[v] = Sherman_MorissonInv(A,u,x)


tmp1 = A\x;
tmp2 = A\u;
v = tmp1 - (tmp2*(u'*tmp1))/(1 + u'*tmp2);


end


% Validation : 
